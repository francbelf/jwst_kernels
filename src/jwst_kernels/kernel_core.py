"""
Uses the Aniano algorithm to generate kernels given a source and target PSF

@author: williams, belfiore
"""

import copy
import warnings

import numpy as np
from astropy.convolution import interpolate_replace_nans, Gaussian2DKernel
from astropy.io import fits
from astropy.modeling import models, fitting
from photutils.centroids import centroid_com
from scipy.ndimage import rotate, zoom
import os

PIXEL_SCALE_NAMES = ['XPIXSIZE', 'CDELT1', 'CD1_1', 'PIXELSCL']

def profile(psf,bins=None, pixscale=1):
    
    i_cen = (psf.shape[0] - 1) / 2
    j_cen = (psf.shape[1] - 1) / 2

    ji, ii = np.meshgrid((np.arange(psf.shape[1]) - j_cen),
                        (np.arange(psf.shape[0]) - i_cen))
    dis = (ji**2 + ii**2)**0.5*pixscale
    
    if bins is None:
        guess_sigma = np.sum(ii[:,int(i_cen)]**2*psf[:,int(i_cen)]/np.sum(psf[:,int(i_cen)]))**0.5*pixscale
        extent = np.min([psf.shape[0]/2*pixscale, guess_sigma*5])
        # int(extent) here used to truncate the bin range to 0 for any
        # sub-arcsec PSF, which made the radial profile (and get_fwhm)
        # silently return zeros. Keep the upper bound in floating point.
        n_bins = max(int(extent/pixscale/2), 8)
        bins = np.linspace(0, extent, n_bins)
            
 
    bin_means = (np.histogram(dis, bins, weights=psf)[0] /
                 np.histogram(dis, bins)[0])

    norm_bins = (bins[:-1]+np.diff(bins)/2)
    return norm_bins, bin_means


def get_fwhm(psf, pixscale=1):

    psf=psf/np.nanmax(psf)
    norm_bins, bin_means = profile(psf,  pixscale=pixscale)

    hwhm = np.interp(0.5, bin_means[::-1]/np.nanmax(bin_means), norm_bins[::-1])
    fwhm = 2*hwhm
    return fwhm

def get_pixscale(hdu):
    """Get pixel scale from header.

    Checks HDU header for pixel scale keywords, and returns a pixel scale in arcsec. If no suitable keyword is found,
    will throw up an error.

    Args:
        hdu (astropy.fits.PrimaryHDU): HDU to get pixel scale for.

    Returns:
        pix_scale (float): Pixel scale in arcsec.

    Raises:
        Warning: If no suitable pixel scale keyword is found in header.

    """
    if isinstance(hdu, (fits.hdu.image.PrimaryHDU, fits.hdu.image.ImageHDU) ):
        header =hdu.header
    elif isinstance(hdu, fits.header.Header ):
        header =hdu
    
    for pixel_keyword in PIXEL_SCALE_NAMES:
        try:
            try:
                pix_scale = np.abs(float(header[pixel_keyword]))
            except ValueError:
                continue
            if pixel_keyword in ['CDELT1', 'CD1_1']:
                pix_scale *= 3600
            return pix_scale
        except KeyError:
            pass

    raise Warning('No pixel scale found')


def fit_2d_gaussian(data, pixscale=None):
    """Fit 2D Gaussian to PSF.

    """

    # Normalize to peak of 1
    data /= np.nanmax(data)

    # Use center of image for first guess center, meshgrid up to feed into modelling
    i_cen = (data.shape[0] - 1) / 2
    j_cen = (data.shape[1] - 1) / 2

    ji, ii = np.meshgrid((np.arange(data.shape[1]) - j_cen),
                         (np.arange(data.shape[0]) - i_cen))
    if pixscale is not None:
        ji *= pixscale
        ii *= pixscale

    guess_sigma = np.sum(ii[:,int(i_cen)]**2*data[:,int(i_cen)]/np.sum(data[:,int(i_cen)]))**0.5
    # Set up model
    model = models.Gaussian2D(amplitude=1,
                              x_mean=0, y_mean=0,
                              x_stddev=guess_sigma, y_stddev=guess_sigma,
                              theta=0
                              )
    fitter = fitting.LevMarLSQFitter()
    fit = fitter(model, ji, ii, data)

    # Take FWHM from mean of the stddevs
    fwhm = 2.355 * (fit.x_stddev.value + fit.y_stddev.value) / 2

    return fwhm


def interp_nans(data, x_stddev=2):
    """Interpolate over any NaNs present in an image.

    Uses astropy.convolution interpolate_replace_nans to smooth over any gaps left in an image.

    Args:
        data (numpy.ndarray): Input data to interpolate NaNs over.
        x_stddev (int, optional): Standard deviation of the Gaussian kernel.
            Defaults to 2 (pixels).

    Returns:
        numpy.ndarray: The data with NaNs interpolated over

    """

    kernel = Gaussian2DKernel(x_stddev=x_stddev)

    image_interp = interpolate_replace_nans(data, kernel)

    return image_interp


def centroid(data):
    i_cen = (data.shape[0] - 1) / 2
    j_cen = (data.shape[1] - 1) / 2

    j_centroid, i_centroid = centroid_com(data)

    # Shift the PSF to center it

    i_shift = int(np.round(i_cen - i_centroid))
    j_shift = int(np.round(j_cen - j_centroid))

    data = np.roll(data, i_shift, axis=0)
    data = np.roll(data, j_shift, axis=1)

    return data


def resample(data, source_pixscale, target_pixscale, interp_order=3):
    """Resample data from one pixel scale to another

    """

    data_resample = zoom(data, source_pixscale / target_pixscale, order=interp_order)

    # force odd-sized array - the kernel needs to be odd

    if data_resample.shape[0] % 2 == 0:
        data_resample = data_resample[:data_resample.shape[0] - 1, :]
    if data_resample.shape[1] % 2 == 0:
        data_resample = data_resample[:, :data_resample.shape[1] - 1]

    return data_resample


def circularize(data, rotations=14):
    """Circularize a PSF.

    Rotate the PSF a number of times, taking an iterative average each time. This serves to make the PSF rotationally
    invariant.

    Args:
        data:
        rotations:

    Returns:

    """

    for n in range(rotations, 0, -1):
        data_rotate = rotate(data, 360 / (2 ** n), order=3, reshape=False)
        data = 0.5 * (data + data_rotate)

    # Set anything outside the maximum radius contained within the whole square to be 0

    radius = np.min(data.shape) / 2
    i_cen = (data.shape[0] - 1) / 2
    j_cen = (data.shape[1] - 1) / 2

    ji, ii = np.meshgrid((np.arange(data.shape[1]) - j_cen),
                         (np.arange(data.shape[0]) - i_cen))

    ri = ji ** 2 + ii ** 2
    data[ri > radius ** 2] = 0

    # Round anything within machine uncertainty to 0

    data[np.abs(data) < np.finfo(float).eps] = 0

    return data


def resize(data, pixscale, grid_size_arcsec=None):
    """Resize data to optimized grid size.

    Args:
        data:
        pixscale:
        grid_size_arcsec:

    Returns:

    TODO:
        * Will fail on non-square arrays

    """

    if grid_size_arcsec is None:
        grid_size_arcsec = np.array([729, 729])

    grid_size_pix = grid_size_arcsec / pixscale
    grid_size_pix = np.array([int(grid_size_pix[0]), int(grid_size_pix[1])])
    #print('inside resize', data.shape, grid_size_pix)
    
    if grid_size_pix[0] % 2 == 0:
        grid_size_pix += 1
    
    if np.all(data.shape > grid_size_pix):
        data_resized = trim(data, grid_size_pix)
    elif np.all(data.shape < grid_size_pix):
        data_resized = zero_pad(data, grid_size_pix)
    elif np.all(data.shape == grid_size_pix):
        data_resized = data
    else:
        raise Warning('Only square arrays implemented!')

    return data_resized


def trim(data, shape):
    """Trim data

    Args:
        data:
        shape:

    Returns:

    """

    shape = np.asarray(shape, dtype=int)
    imshape = np.asarray(data.shape, dtype=int)
    dshape = imshape - shape

    if np.all(imshape == shape):
        return data

    if np.any(dshape % 2 != 0):
        raise ValueError('Source and target shapes have different parity.')

    idx, idy = np.indices(shape)
    offx, offy = dshape // 2

    return data[idx + offx, idy + offy]


def zero_pad(data, shape):
    """Zero pad data

    Args:
        data:
        shape:

    Returns:

    """

    shape = np.asarray(shape, dtype=int)
    imshape = np.asarray(data.shape, dtype=int)
    dshape = shape - imshape

    if np.all(imshape == shape):
        return data

    data_pad = np.zeros(shape)
    idx, idy = np.indices(imshape)

    if np.any(dshape % 2 != 0):
        raise ValueError('Source and target shapes have different parity.')
    offx, offy = dshape // 2

    data_pad[idx + offx, idy + offy] = data

    return data_pad


def high_pass_filter(data, fwhm, pixscale=0.1):
    """Short desc

    Args:
        data (np.ndarray): Data to high-pass filter.
        fwhm (float): FWHM, in arcsec.
        pixscale (float, optional): Pixel scale, in arcsec. Defaults to 0.1.

    Returns:

    """

    k_b = 8 * np.pi / (fwhm * pixscale)
    k_a = 0.9 * k_b
    i_cen = (data.shape[0] - 1) / 2
    j_cen = (data.shape[1] - 1) / 2

    ji, ii = np.meshgrid((np.arange(data.shape[1]) - j_cen),
                         (np.arange(data.shape[0]) - i_cen))
    ri = np.sqrt(ji ** 2 + ii ** 2) * pixscale

    # Create filter
    data_filter = np.zeros(data.shape)
    data_filter[ri <= k_a] = 1

    idx = np.where((k_a < ri) & (ri <= k_b))
    data_filter[idx] = np.exp(-(1.8249 * (ri[idx] - k_a) / (k_b - k_a)) ** 4)

    # Apply filter
    data_filtered = data_filter * data

    return data_filtered


def low_pass_filter(data, pixscale=0.1):
    """Low-pass filter.

    """
    # Calculate k_h as the first value for which the radial value of the FFT drops below 0.5% of the max value

    i_range = data.shape[0]
    j_range = data.shape[1]

    data_slice = data[
                 int((i_range - 1) / 2):,
                 int((j_range - 1) / 2)
                 ]
    data_slice_max = np.nanmax(data_slice)

    try:
        k_h = np.where(data_slice < 0.005 * data_slice_max)[0][0] * pixscale
    except IndexError:
        print('k_h too large. Something is probably up with kernel generation...')
        k_h = len(data_slice) * pixscale
    k_l = 0.7 * k_h

    i_cen = (data.shape[0] - 1) / 2
    j_cen = (data.shape[1] - 1) / 2

    ji, ii = np.meshgrid((np.arange(data.shape[1]) - j_cen),
                         (np.arange(data.shape[0]) - i_cen))
    ri = np.sqrt(ji ** 2 + ii ** 2) * pixscale

    # Create the low-pass filter
    data_filter = np.zeros(data.shape)

    idx = np.where((k_l < ri) & (ri <= k_h))
    data_filter[idx] = 0.5 * (1 + np.cos(np.pi * (ri[idx] - k_l) / (k_h - k_l)))

    data_filter[ri <= k_l] = 1

    return data_filter


def trim_kernel_energy(kernel, energy_tol=0.999):
    """ Trim kernel based on enclosed energy to speed up later convolutions/space requirements

    Args:
        kernel:
        energy_tol:

    Returns:

    TODO:
        * This could fail if the kernel is too small/too big. Keep an eye out

    """
    kernel_radius = int((kernel.shape[0] - 1) / 2)

    i_cen = (kernel.shape[0] - 1) / 2
    j_cen = (kernel.shape[1] - 1) / 2
    ji, ii = np.meshgrid((np.arange(kernel.shape[1]) - j_cen),
                         (np.arange(kernel.shape[0]) - i_cen))

    ri = np.sqrt(ji ** 2 + ii ** 2)

    total_kernel_energy = np.nansum(np.abs(kernel[ri <= kernel_radius]))

    for radius in range(kernel_radius):
        idx = np.where(ri <= radius)
        enclosed_energy = np.nansum(np.abs(kernel[idx]))
        frac_kernel_energy = enclosed_energy / total_kernel_energy
        if frac_kernel_energy >= energy_tol:
            break
    trim_shape = (radius * 2 + 1, radius * 2 + 1)

    kernel_trimmed = trim(kernel, trim_shape)

    return kernel_trimmed


class MakeConvolutionKernel:
    """Class to generate kernels following the Aniano 2011 algorithm.
    """

    def __init__(self,
                 source_psf=None,
                 source_fwhm=None,
                 source_name=None,
                 source_pixscale=1,
                 target_psf=None,
                 target_fwhm=None,
                 target_name=None,
                 target_pixscale=1,
                 common_pixscale=None,
                 grid_size_arcsec=None,
                 verbose=False,
                 ):
        """
        Args:
            source_psf: Source PSF as ndarray or FITS HDU. Required.
            source_fwhm: Source FWHM in arcsec. Fitted if not provided.
            source_name: Label for source PSF. Defaults to 'source'.
            source_pixscale: Pixel scale in arcsec (used when source_psf is an array).
            target_psf: Target PSF as ndarray or FITS HDU. Optional -- not needed
                if only processing the source PSF via process_source_psf().
            target_fwhm: Target FWHM in arcsec. Fitted if not provided.
            target_name: Label for target PSF. Defaults to 'target'.
            target_pixscale: Pixel scale in arcsec (used when target_psf is an array).
            common_pixscale: Common pixel scale to resample both PSFs onto.
                Defaults to source_pixscale if not provided.
            grid_size_arcsec: Grid size in arcsec for the resize step.
                Defaults to the source PSF's native grid size if not provided.
            verbose: Print progress messages.
        """
        if source_psf is None:
            raise ValueError('source_psf must be provided')

        if isinstance(source_psf, (fits.hdu.image.PrimaryHDU, fits.hdu.image.ImageHDU)):
            self.source_psf = copy.deepcopy(source_psf.data)
            self.source_pixscale = get_pixscale(source_psf)
        elif isinstance(source_psf, np.ndarray):
            self.source_psf = copy.deepcopy(source_psf)
            self.source_pixscale = source_pixscale
        else:
            raise TypeError('source_psf must be an ndarray or FITS HDU')

        self.target_psf = None
        self.target_pixscale = target_pixscale
        self.target_fwhm = None

        if target_psf is not None:
            if isinstance(target_psf, (fits.hdu.image.PrimaryHDU, fits.hdu.image.ImageHDU)):
                self.target_psf = copy.deepcopy(target_psf.data)
                self.target_pixscale = get_pixscale(target_psf)
            elif isinstance(target_psf, np.ndarray):
                self.target_psf = copy.deepcopy(target_psf)
                self.target_pixscale = target_pixscale
            else:
                raise TypeError('target_psf must be an ndarray or FITS HDU')

            if not target_fwhm:
                target_fwhm = fit_2d_gaussian(data=self.target_psf, pixscale=self.target_pixscale)
            self.target_fwhm = copy.deepcopy(target_fwhm)

        if not source_fwhm:
            source_fwhm = fit_2d_gaussian(data=self.source_psf, pixscale=self.source_pixscale)

        if self.target_fwhm is not None and source_fwhm >= self.target_fwhm:
            raise Warning('Cannot create kernel from lower to higher resolution data!', source_fwhm, self.target_fwhm)

        self.source_fwhm = copy.deepcopy(source_fwhm)

        self.source_name = copy.deepcopy(source_name) if source_name is not None else 'source'
        self.target_name = copy.deepcopy(target_name) if target_name is not None else 'target'


        self.common_pixscale = copy.deepcopy(common_pixscale)
        self.grid_size_arcsec = copy.deepcopy(grid_size_arcsec) if grid_size_arcsec is not None else None

        self.source_fourier = None
        self.target_fourier = None

        self.kernel_fourier = None
        self.kernel = None

        # Flags to track if PSFs have been processed (in case process_source_psf or make_convolution_kernel are called multiple times)
        self._source_processed = False
        self._target_processed = False

        self.verbose = verbose

    def _process_psf(self, psf_data, pixscale, common_pixscale, grid_size_arcsec,
                     do_circularize=True):
        """Spatially process a single PSF: interp NaNs, resample, centroid, circularize, resize, normalize.

        If ``do_circularize`` is False, the real-space rotate-and-average step is skipped.
        """
        psf_data = interp_nans(psf_data)

        if not np.isclose(pixscale, common_pixscale):
            if self.verbose:
                print('Resampling')
            psf_data = resample(psf_data, pixscale, common_pixscale)
            if self.verbose:
                print(f'  resampled shape: {psf_data.shape}')
        else:
            if self.verbose:
                print('Skipping resample (pixscales match)')
            # resample() forces odd dimensions; do the same when skipping
            if psf_data.shape[0] % 2 == 0:
                psf_data = psf_data[:psf_data.shape[0] - 1, :]
            if psf_data.shape[1] % 2 == 0:
                psf_data = psf_data[:, :psf_data.shape[1] - 1]

        if self.verbose:
            print('Centroiding')
        psf_data = centroid(psf_data)

        if do_circularize:
            if self.verbose:
                print('Circularizing')
            psf_data = circularize(psf_data)
        elif self.verbose:
            print('Skipping real-space circularize (do_circularize=False)')

        if self.verbose:
            print('Resizing')
            if grid_size_arcsec is not None:
                print('  grid in pixels', grid_size_arcsec / common_pixscale)
        psf_data = resize(psf_data, common_pixscale, grid_size_arcsec=grid_size_arcsec)

        psf_data /= np.nansum(psf_data)

        if self.verbose:
            print(f'  final shape: {psf_data.shape}')

        return psf_data

    def process_source_psf(self, do_circularize=True, do_fourier_filter=True):
        """Spatially process the source PSF and optionally apply Fourier-domain filtering.

        Applies the Aniano 2011 algorithm to the source PSF: interp NaNs,
        resample, centroid, (optionally) real-space circularize, resize,
        normalize, then (optionally) FFT -> circularize -> high-pass
        filter -> IFFT.

        Args:
            do_circularize: If True (default), apply the real-space rotate-and-average
                step during spatial processing. If False, skip it.
            do_fourier_filter: If True (default), apply the Fourier-domain
                circularize + high-pass filter block. If False, return after
                spatial processing only.

        Variants:
            * ``do_circularize=True,  do_fourier_filter=True``  -> full Aniano processing.
            * ``do_circularize=True,  do_fourier_filter=False`` -> real-space circularized only.
            * ``do_circularize=False, do_fourier_filter=True``  -> Fourier-filtered without
              real-space circularization.
            * ``do_circularize=False, do_fourier_filter=False`` -> spatial processing only
              (interp NaNs, resample, centroid, resize, normalize). Neither the real-space
              rotate-and-average circularize nor the Fourier-domain circularize + high-pass
              filter block is applied. 

        If common_pixscale is None, defaults to source_pixscale (no resample).
        If grid_size_arcsec is None, defaults to the source PSF's native grid
        size so that no resize is needed either.

        After calling, the processed source PSF is available via
        self.source_psf and can be saved with save_processed_psf().

        """
        if self._source_processed:
            raise RuntimeError('Source PSF has already been processed. '
                               'Create a new MakeConvolutionKernel instance to reprocess.')

        if self.common_pixscale is None:
            self.common_pixscale = self.source_pixscale
        if self.grid_size_arcsec is None:
            self.grid_size_arcsec = (
                np.array(self.source_psf.shape, dtype=float) * self.source_pixscale
            )

        if self.verbose:
            print(f'Processing source PSF ({self.source_name})')

        self.source_psf = self._process_psf(
            self.source_psf, self.source_pixscale,
            self.common_pixscale, self.grid_size_arcsec,
            do_circularize=do_circularize,
        )

        if do_fourier_filter:
            if self.verbose:
                print('Fourier-domain filtering')

            source_fourier = np.real(np.fft.fft2(np.fft.ifftshift(self.source_psf)))
            source_fourier = np.fft.fftshift(source_fourier)
            source_fourier = circularize(source_fourier)
            source_fourier = high_pass_filter(
                source_fourier, self.source_fwhm, self.common_pixscale / 2,
            )
            self.source_psf = np.fft.fftshift(
                np.real(np.fft.ifft2(np.fft.ifftshift(source_fourier)))
            )
            self.source_psf[np.abs(self.source_psf) <= np.finfo(float).eps] = 0
            self.source_psf /= np.nansum(self.source_psf)
        elif self.verbose:
            print('Skipping Fourier-domain filtering (do_fourier_filter=False)')

        self._source_processed = True

    def make_convolution_kernel(self):
        """Generate a convolution kernel from source to target PSF using the Aniano algorithm.

        Requires that target_psf was provided at init time.

        If common_pixscale is None, defaults to source_pixscale.
        If grid_size_arcsec is None, resize() defaults to [729, 729] arcsec.
        """
        if self.target_psf is None:
            raise ValueError('target_psf is required for make_convolution_kernel. '
                             'Use process_source_psf() to process the source PSF alone.')
        if self._source_processed or self._target_processed:
            raise RuntimeError('PSFs have already been processed. '
                               'Create a new MakeConvolutionKernel instance to reprocess.')

        if self.common_pixscale is None:
            self.common_pixscale = self.source_pixscale

        if self.verbose:
            print('Processing source and target PSFs')

        self.source_psf = self._process_psf(
            self.source_psf, self.source_pixscale,
            self.common_pixscale, self.grid_size_arcsec,
        )
        self.target_psf = self._process_psf(
            self.target_psf, self.target_pixscale,
            self.common_pixscale, self.grid_size_arcsec,
        )
        self._source_processed = True
        self._target_processed = True

        if self.verbose:
            print('common pixel grid', self.source_psf.shape)

        self._fourier_kernel_pipeline()

    def _fourier_kernel_pipeline(self):
        """Run the Fourier-domain kernel construction on already-spatially-processed PSFs.

        Assumes ``self.source_psf`` and ``self.target_psf`` are on the common
        pixel grid, centroided, normalized, and ready for FFT. Populates
        ``self.source_fourier``, ``self.target_fourier``, ``self.kernel_fourier``,
        and ``self.kernel``.
        """
        if self.target_psf is None:
            raise ValueError('target_psf is required to build a convolution kernel.')
        if self.target_fwhm is None:
            raise ValueError('target_fwhm is required to build a convolution kernel.')

        if self.verbose:
            print('FFTing')

        # We now move onto the FFT part. Fourier transform the PSFs - only take the real part

        self.source_fourier = np.real(np.fft.fft2(np.fft.ifftshift(self.source_psf)))
        self.target_fourier = np.real(np.fft.fft2(np.fft.ifftshift(self.target_psf)))

        # Make sure center of FFT is in middle

        self.source_fourier = np.fft.fftshift(self.source_fourier)
        self.target_fourier = np.fft.fftshift(self.target_fourier)

        if self.verbose:
            print('Circularizing')

        # Circularize
        self.source_fourier = circularize(self.source_fourier)
        self.target_fourier = circularize(self.target_fourier)

        if self.verbose:
            print('High-pass filter')

        # High-pass filter
        source_fourier_high_pass = high_pass_filter(self.source_fourier, self.source_fwhm, self.common_pixscale / 2)
        target_fourier_high_pass = high_pass_filter(self.target_fourier, self.target_fwhm, self.common_pixscale / 2)

        if self.verbose:
            print('Inverting')

        # Invert the source fourier, any infs go to 0
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            source_fourier_invert = source_fourier_high_pass ** -1
            inf_idx = np.where(~np.isfinite(source_fourier_invert))
            source_fourier_invert[inf_idx] = 0

        if self.verbose:
            print('Low-pass filter')

        # Low-pass filter
        source_fourier_low_pass = low_pass_filter(self.source_fourier, self.common_pixscale / 2)

        if self.verbose:
            print('Creating kernel')

        # FFT of convolution kernel
        self.kernel_fourier = target_fourier_high_pass * (source_fourier_low_pass * source_fourier_invert)
        self.kernel_fourier = np.fft.ifftshift(self.kernel_fourier)

        if self.verbose:
            print('IFFT-ing')

        # IFFT to kernel and round out any tiny computational errors
        self.kernel = np.fft.fftshift(np.real(np.fft.ifft2(self.kernel_fourier)))
        self.kernel[np.abs(self.kernel) <= np.finfo(float).eps] = 0

        if self.verbose:
            print('Centroiding')
            print('size kernel', self.kernel.shape)

        # Centroid again, just in case
        self.kernel = centroid(self.kernel)

        if self.verbose:
            print('Trimming kernel')

        # Trim kernel based on enclosed energy
        self.kernel = trim_kernel_energy(self.kernel)

        if self.verbose:
            print('Last little bits')

        # Finally, circularize kernel and normalize to peak of 1
        self.kernel = circularize(self.kernel)
        self.kernel /= np.nanmax(self.kernel)

    def make_convolution_kernel_from_processed(self):
        """Generate a convolution kernel from already-processed source and target PSFs.

        Both ``self.source_psf`` and ``self.target_psf`` must already be
        spatially processed (interp NaNs, resample, centroid, optionally
        circularize, resize, normalize) and live on a common pixel grid.
        This method skips the spatial ``_process_psf`` step and runs only
        the Fourier-domain kernel construction pipeline.

        Requires:
            - ``self.target_psf`` is not None.
            - ``self.source_psf.shape == self.target_psf.shape``.
            - ``self.source_pixscale`` is close to ``self.target_pixscale``.
        """
        if self.target_psf is None:
            raise ValueError('target_psf is required for make_convolution_kernel_from_processed. '
                             'Use process_source_psf() to process the source PSF alone.')
        if self._source_processed or self._target_processed:
            raise RuntimeError('PSFs have already been processed. '
                               'Create a new MakeConvolutionKernel instance to reprocess.')
        if self.source_psf.shape != self.target_psf.shape:
            raise ValueError(
                f'source_psf and target_psf must have the same shape, got '
                f'{self.source_psf.shape} vs {self.target_psf.shape}.'
            )
        if not np.isclose(self.source_pixscale, self.target_pixscale):
            raise ValueError(
                f'source_pixscale ({self.source_pixscale}) and target_pixscale '
                f'({self.target_pixscale}) must be (close to) equal.'
            )

        if self.common_pixscale is None:
            self.common_pixscale = self.source_pixscale
        elif not np.isclose(self.common_pixscale, self.source_pixscale):
            raise ValueError(
                f'common_pixscale ({self.common_pixscale}) must match '
                f'source/target pixscale ({self.source_pixscale}).'
            )

        self._source_processed = True
        self._target_processed = True

        if self.verbose:
            print('Skipping spatial processing (inputs already processed)')
            print('common pixel grid', self.source_psf.shape)

        self._fourier_kernel_pipeline()

    def save_processed_psf(self, outdir, which='source', filename_suffix='aniano_circ_filt'):
        """Save a processed PSF to a FITS file.

        Args:
            outdir (str): Output directory.
            which (str): Which PSF to save. One of 'source', 'target', or 'both'.
                Defaults to 'source'.

        Returns:
            list[str]: Paths to the saved FITS files.
        """
        if which not in ('source', 'target', 'both'):
            raise ValueError(f"which must be 'source', 'target', or 'both', got '{which}'")

        os.makedirs(outdir, exist_ok=True)
        saved = []

        psfs_to_save = []
        if which in ('source', 'both'):
            psfs_to_save.append((self.source_psf, self.source_name))
        if which in ('target', 'both'):
            psfs_to_save.append((self.target_psf, self.target_name))

        for psf_data, name in psfs_to_save:
            if psf_data is None:
                raise ValueError(f"No PSF data for '{name}'. Process the PSF first.")

            hdu = fits.PrimaryHDU(data=np.array(psf_data, dtype=np.float32))
            hdu.header['PSFNAME'] = (name, 'PSF identifier')
            hdu.header['PIXELSCL'] = (self.common_pixscale, 'arcsec/pixel')
            hdu.header['CRPIX1'] = (psf_data.shape[1] + 1) / 2
            hdu.header['CRPIX2'] = (psf_data.shape[0] + 1) / 2
            hdu.header['CRVAL1'] = 0.0
            hdu.header['CRVAL2'] = 0.0
            hdu.header['CDELT1'] = -self.common_pixscale / 3600
            hdu.header['CDELT2'] = self.common_pixscale / 3600

            file_name = os.path.join(outdir, f'{name}_{filename_suffix}.fits')
            hdu.writeto(file_name, overwrite=True)
            saved.append(file_name)

        return saved

    def write_out_kernel(self, outdir=None, add_keys=None, naming_convention='PHANGS', print_name=False):
        """

        Returns:

        """
        if naming_convention=='PHANGS':
            target_name = self.target_name.replace('.', 'p')
            target_name = target_name.lower()
            source_name = self.source_name.lower()
            file_name = os.path.join( outdir, '%s_to_%s.fits' % (source_name, target_name))
        else:
            file_name = os.path.join( outdir, '%s_to_%s.fits' % (self.source_name, self.target_name))

        # Build the fits file. Use 32bit precision to cut down space

        hdu = fits.PrimaryHDU(data=np.array(self.kernel, dtype=np.float32))

        hdu.header['BITPIX'] = -32

        hdu.header['CRPIX1'] = (self.kernel.shape[1] + 1) / 2
        hdu.header['CRPIX2'] = (self.kernel.shape[0] + 1) / 2

        hdu.header['CRVAL1'] = 0.00
        hdu.header['CRVAL2'] = 0.00

        hdu.header['CDELT1'] = - self.common_pixscale / 3600
        hdu.header['CDELT2'] = self.common_pixscale / 3600

        hdu.header['CTYPE1'] = 'RA---TAN'
        hdu.header['CTYPE2'] = 'DEC--TAN'

        if (add_keys is not None) & (type(add_keys)==dict):
            for key in  add_keys:
                hdu.header[key] = add_keys[key]

        hdu.writeto(file_name, overwrite=True)
        if print_name:
            print(f"Saved: {file_name}")



def make_jwst_cross_kernel(input_filter, target_filter, psf_dir=None,
                           outdir=None, save_kernel=False,
                           common_pixscale=None, detector_effects=True,
                           naming_convention='PHANGS', verbose=False,
                           save_plots=False, plot_dir=None):
    """Generate the kernel that matches a JWST input filter to a target JWST filter.

    Works for both MIRI and NIRCam. If the source or target raw PSF is not
    found under ``psf_dir``, ``read_PSF`` will regenerate it via webbpsf.

    Parameters
    ----------
    input_filter, target_filter : dict
        Each must have ``'camera'`` and ``'filter'`` keys.
    psf_dir, outdir : str, optional
        Where raw PSFs live, and where kernels are written.
    save_kernel : bool
        Whether to write the kernel to disk.
    common_pixscale : float, optional
        Force the common pixel scale used for resampling. Defaults to the
        source pixel scale.
    detector_effects : bool
        Whether to pull the ``OVERDIST`` extension (detector effects on).
    naming_convention : str
        Filename convention passed through to ``write_out_kernel``.
    save_plots : bool
        If True, also write a diagnostic PNG via
        :func:`jwst_kernels.evaluate_kernels.plot_kernel_diagnostic`.
    plot_dir : str, optional
        Directory for diagnostic PNGs. Defaults to ``<outdir>/plots`` when
        ``save_plots=True`` and ``plot_dir`` is None.

    Returns
    -------
    MakeConvolutionKernel
    """
    from jwst_kernels.make_psf import read_PSF

    source_psf, source_pixscale = read_PSF(
        input_filter, detector_effects=detector_effects, psf_dir=psf_dir)
    target_psf, target_pixscale = read_PSF(
        target_filter, detector_effects=detector_effects, psf_dir=psf_dir)

    if common_pixscale is None:
        common_pixscale = source_pixscale

    grid_size_arcsec = np.array([361 * common_pixscale,
                                 361 * common_pixscale])

    kk = MakeConvolutionKernel(source_psf=source_psf,
                               source_pixscale=source_pixscale,
                               source_name=input_filter['filter'],
                               target_psf=target_psf,
                               target_pixscale=target_pixscale,
                               target_name=target_filter['filter'],
                               common_pixscale=common_pixscale,
                               grid_size_arcsec=grid_size_arcsec,
                               verbose=verbose)
    kk.make_convolution_kernel()
    dict_extension = {'DETEF': detector_effects}
    if outdir is None:
        outdir = os.path.abspath('.')
    if save_kernel:
        kk.write_out_kernel(outdir=outdir, add_keys=dict_extension,
                            naming_convention=naming_convention)

    if save_plots:
        from jwst_kernels.evaluate_kernels import plot_kernel_diagnostic
        if plot_dir is None:
            plot_dir = os.path.join(str(outdir), 'plots')
        plot_kernel_diagnostic(kk, plot_dir=plot_dir)

    return kk


def make_jwst_kernel_to_Gauss(input_filter, target_gaussian, psf_dir=None,
                              outdir=None, save_kernel=False,
                              detector_effects=True,
                              naming_convention='PHANGS', verbose=False,
                              size_kernel_asec=None,
                              save_plots=False, plot_dir=None):
    """Generate the kernel that matches a JWST filter PSF to a Gaussian target.

    Parameters
    ----------
    input_filter : dict
        Must have ``'camera'`` and ``'filter'`` keys.
    target_gaussian : dict
        Must have a ``'fwhm'`` key (arcsec).
    psf_dir, outdir : str, optional
    save_kernel : bool
    detector_effects : bool
    naming_convention : str
    size_kernel_asec : float, optional
        Override the auto-chosen kernel footprint (arcsec).
    save_plots : bool
        If True, also write a diagnostic PNG via
        :func:`jwst_kernels.evaluate_kernels.plot_kernel_diagnostic`.
    plot_dir : str, optional
        Directory for diagnostic PNGs. Defaults to ``<outdir>/plots`` when
        ``save_plots=True`` and ``plot_dir`` is None.

    Returns
    -------
    MakeConvolutionKernel
    """
    from jwst_kernels.make_psf import read_PSF, makeGaussian_2D

    source_psf, source_pixscale = read_PSF(
        input_filter, detector_effects=detector_effects, psf_dir=psf_dir)

    common_pixscale = source_pixscale
    target_pixscale = source_pixscale

    if size_kernel_asec is None:
        sz = int(target_gaussian['fwhm'] / 2.355 * 20 / target_pixscale)
    else:
        sz = int(size_kernel_asec / target_pixscale)
    if sz % 2 == 0:
        sz += 1

    yy, xx = np.meshgrid(np.arange(sz) - (sz - 1) / 2,
                         np.arange(sz) - (sz - 1) / 2)
    target_psf = makeGaussian_2D(
        (xx, yy), (0, 0),
        (target_gaussian['fwhm'] / 2.355 / target_pixscale,
         target_gaussian['fwhm'] / 2.355 / target_pixscale))
    target_name = 'gauss{:.2f}'.format(target_gaussian['fwhm'])

    if verbose:
        print('making kernel from ' + input_filter['filter']
              + ' with source PSF to Gaussian with FWHM '
              + '{:.3f}'.format(target_gaussian['fwhm']))

    grid_size_arcsec = np.array([sz * target_pixscale, sz * target_pixscale])

    kk = MakeConvolutionKernel(source_psf=source_psf,
                               source_pixscale=source_pixscale,
                               source_name=input_filter['filter'],
                               target_psf=target_psf,
                               target_pixscale=target_pixscale,
                               target_name=target_name,
                               common_pixscale=common_pixscale,
                               grid_size_arcsec=grid_size_arcsec,
                               verbose=verbose)
    kk.make_convolution_kernel()
    dict_extension = {'DETEF': detector_effects}
    if outdir is None:
        outdir = os.path.abspath('.')
    if save_kernel:
        kk.write_out_kernel(outdir=outdir, add_keys=dict_extension,
                            naming_convention=naming_convention)

    if save_plots:
        from jwst_kernels.evaluate_kernels import plot_kernel_diagnostic
        if plot_dir is None:
            plot_dir = os.path.join(str(outdir), 'plots')
        plot_kernel_diagnostic(kk, plot_dir=plot_dir)

    return kk


def plot_kernel(kk, save_plot=False, save_dir=None, want_convolve=True):
    """Plot source PSF, target PSF, and the kernel with radial profiles.

    Parameters
    ----------
    kk : MakeConvolutionKernel
    save_plot : bool
    save_dir : str, optional
    want_convolve : bool
        Also convolve source x kernel and overlay it on the radial profile
        (can be slow for large kernels).
    """
    import matplotlib.pyplot as plt
    from astropy.convolution import convolve

    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(12, 4))

    ax1.imshow(np.log10(kk.source_psf / np.max(kk.source_psf)),
               vmax=0, vmin=-4)
    ax1.set_title(kk.source_name)

    ax2.imshow(np.log10(kk.target_psf / np.max(kk.target_psf)),
               vmax=0, vmin=-4)
    ax2.set_title(kk.target_name)

    extent = int(10 * kk.target_fwhm / kk.common_pixscale / 2)
    ax3.plot(*profile(kk.source_psf / np.max(kk.source_psf),
                      bins=np.linspace(0, 6 * kk.target_fwhm, extent),
                      pixscale=kk.common_pixscale),
             c='b', label=kk.source_name)
    ax3.plot(*profile(kk.target_psf / np.max(kk.target_psf),
                      bins=np.linspace(0, 6 * kk.target_fwhm, extent),
                      pixscale=kk.common_pixscale),
             c='k', label=kk.target_name, lw=5)
    ax3.plot(*profile(kk.kernel / np.max(kk.kernel),
                      bins=np.linspace(0, 6 * kk.target_fwhm, extent),
                      pixscale=kk.common_pixscale),
             c='g', label='kernel')
    if want_convolve:
        target_conv = convolve(kk.source_psf, kk.kernel)
        ax3.plot(*profile(target_conv / np.max(target_conv),
                          bins=np.linspace(0, 6 * kk.target_fwhm, extent),
                          pixscale=kk.common_pixscale),
                 c='r', label='model', ls='-')

        xx, targetprof = profile(
            kk.target_psf / np.max(kk.target_psf),
            bins=np.linspace(0, 6 * kk.target_fwhm, extent),
            pixscale=kk.common_pixscale)
        _, convprof = profile(
            target_conv / np.max(target_conv),
            bins=np.linspace(0, 6 * kk.target_fwhm, extent),
            pixscale=kk.common_pixscale)
        ax3.plot(xx, (targetprof - convprof),
                 c='r', label='residual', ls='--')
    ax3.legend()
    ax3.set_ylim([-0.1, 1.1])
    ax3.set_xlim([0, 6 * kk.target_fwhm])

    if save_plot:
        if save_dir is None:
            save_dir = os.path.abspath('.')
        plt.savefig(os.path.join(save_dir,
                                 kk.source_name + '_' + kk.target_name + '.png'),
                    dpi=300)
