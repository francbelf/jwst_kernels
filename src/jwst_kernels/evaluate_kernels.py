#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits, ascii
import astropy.table as table
from jwst_kernels.kernel_core import (
    fit_2d_gaussian,
    get_fwhm,
    make_jwst_kernel_to_Gauss,
    plot_kernel,
    profile,
)
from jwst_kernels.make_psf import  read_PSF
from astropy.convolution import convolve, convolve_fft
from os import path
import jwst_kernels
from scipy import interpolate

def evaluate_kernel(kk):
    
    kk.kernel=kk.kernel/np.sum(kk.kernel)
    target_conv = convolve(kk.source_psf, kk.kernel)
    # D kernel performance measure Aniano Eq 20
    D = np.sum(np.abs(target_conv-kk.target_psf))
    # Wm kernel performance measure Aniano eq 21
    Wm = 0.5 *np.sum( np.abs(kk.kernel) - kk.kernel)
    return D, Wm


def plot_kernel_diagnostic(kk, plot_dir=None, save_path=None,
                           want_convolve=True, dpi=150):
    """Save a 2x2 diagnostic plot for a constructed convolution kernel.

    Layout
    ------
    (0, 0) Source PSF, log10 scale, vmin=-4.
    (0, 1) Target PSF, log10 scale, vmin=-4.
    (1, 0) Kernel image with a signed colormap so positive and negative
           lobes are both visible.
    (1, 1) Radial profiles of the source PSF, target PSF, and
           ``source * kernel`` (the matched model), with a stats text box
           reporting the Aniano D (Eq. 20) and W_- (Eq. 21) performance
           measures together with source / target / recovered FWHMs.

    Calls :func:`evaluate_kernel`, which renormalises ``kk.kernel`` to
    ``sum=1`` in place.

    Parameters
    ----------
    kk : MakeConvolutionKernel
        Kernel object with ``kernel``, ``source_psf``, and ``target_psf``
        populated.
    plot_dir : str, optional
        Directory to save the PNG into. The filename is derived from
        ``kk.source_name`` / ``kk.target_name`` using the PHANGS lowercase
        ``"."`` -> ``"p"`` convention. Ignored if ``save_path`` is given.
    save_path : str, optional
        Explicit output PNG path. Overrides ``plot_dir``.
    want_convolve : bool
        Also compute ``source * kernel`` and overlay its radial profile
        plus a residual curve (target - model).
    dpi : int
        Figure DPI for saved PNGs.

    Returns
    -------
    str or matplotlib.figure.Figure
        Path to the saved PNG when a save destination was given (figure is
        closed). Otherwise the open figure for the caller to display.
    """
    if kk.kernel is None:
        raise ValueError("kk.kernel is None; build the kernel before plotting.")
    if kk.source_psf is None or kk.target_psf is None:
        raise ValueError("source_psf/target_psf are required for plotting.")

    # Compute Aniano D (Eq. 20) and W_- (Eq. 21) inline. Using FFT
    # convolution avoids the very slow spatial path on ~361x361 arrays and
    # lets us reuse target_conv for the radial profile in panel (1, 1).
    kk.kernel = kk.kernel / np.sum(kk.kernel)
    target_conv = convolve_fft(kk.source_psf, kk.kernel,
                               normalize_kernel=False,
                               allow_huge=True)
    D = float(np.sum(np.abs(target_conv - kk.target_psf)))
    Wm = float(0.5 * np.sum(np.abs(kk.kernel) - kk.kernel))

    source_fwhm = get_fwhm(kk.source_psf, pixscale=kk.common_pixscale)
    if kk.target_fwhm is not None:
        target_fwhm = float(kk.target_fwhm)
    else:
        target_fwhm = get_fwhm(kk.target_psf, pixscale=kk.common_pixscale)

    target_conv_fwhm = None
    if want_convolve:
        target_conv_fwhm = get_fwhm(target_conv, pixscale=kk.common_pixscale)

    fig, axes = plt.subplots(2, 2, figsize=(10, 9))

    src_norm = kk.source_psf / np.nanmax(kk.source_psf)
    axes[0, 0].imshow(np.log10(np.clip(src_norm, 1e-6, None)),
                      vmin=-4, vmax=0, origin='lower')
    axes[0, 0].set_title(f"{kk.source_name}\n(FWHM={source_fwhm:.3f}\")")
    axes[0, 0].set_xlabel('pixel')
    axes[0, 0].set_ylabel('pixel')

    tgt_norm = kk.target_psf / np.nanmax(kk.target_psf)
    axes[0, 1].imshow(np.log10(np.clip(tgt_norm, 1e-6, None)),
                      vmin=-4, vmax=0, origin='lower')
    axes[0, 1].set_title(f"{kk.target_name}\n(FWHM={target_fwhm:.3f}\")")
    axes[0, 1].set_xlabel('pixel')
    axes[0, 1].set_ylabel('pixel')

    kabs = np.nanmax(np.abs(kk.kernel))
    if kabs == 0 or not np.isfinite(kabs):
        kabs = 1.0
    knorm = kk.kernel / kabs
    axes[1, 0].imshow(knorm, cmap='RdBu_r', vmin=-1, vmax=1, origin='lower')
    axes[1, 0].set_title('kernel\n(red=positive, blue=negative)')
    axes[1, 0].set_xlabel('pixel')
    axes[1, 0].set_ylabel('pixel')

    extent = max(int(10 * target_fwhm / kk.common_pixscale / 2), 20)
    bins = np.linspace(0, 6 * target_fwhm, extent)

    src_x, src_y = profile(kk.source_psf / np.nanmax(kk.source_psf),
                           bins=bins, pixscale=kk.common_pixscale)
    tgt_x, tgt_y = profile(kk.target_psf / np.nanmax(kk.target_psf),
                           bins=bins, pixscale=kk.common_pixscale)

    ax_p = axes[1, 1]
    ax_p.plot(src_x, src_y, c='b', label=kk.source_name)
    ax_p.plot(tgt_x, tgt_y, c='k', label=kk.target_name, lw=3)

    if want_convolve and target_conv is not None:
        mod_x, mod_y = profile(target_conv / np.nanmax(target_conv),
                               bins=bins, pixscale=kk.common_pixscale)
        ax_p.plot(mod_x, mod_y, c='r', label='source * kernel')
        ax_p.plot(mod_x, tgt_y - mod_y, c='r', ls='--',
                  label='target - model')

    ax_p.set_xlabel('radius (arcsec)')
    ax_p.set_ylabel('normalised intensity')
    ax_p.set_xlim([0, 6 * target_fwhm])
    ax_p.set_ylim([-0.1, 1.1])
    ax_p.legend(loc='upper right', fontsize=8)

    stats_lines = [
        f"D = {D:.3e}",
        f"W$_-$ = {Wm:.3e}",
        f"FWHM source = {source_fwhm:.3f}\"",
        f"FWHM target = {target_fwhm:.3f}\"",
    ]
    if target_conv_fwhm is not None:
        stats_lines.append(f"FWHM source*kernel = {target_conv_fwhm:.3f}\"")

    ax_p.text(0.02, 0.98, "\n".join(stats_lines),
              transform=ax_p.transAxes,
              ha='left', va='top', fontsize=9,
              bbox=dict(boxstyle='round', facecolor='white', alpha=0.85,
                        edgecolor='gray'))

    fig.suptitle(f"{kk.source_name} -> {kk.target_name}", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    out_path = None
    if save_path is not None:
        out_path = str(save_path)
        os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
    elif plot_dir is not None:
        target_name = kk.target_name.replace('.', 'p').lower()
        source_name = kk.source_name.lower()
        os.makedirs(str(plot_dir), exist_ok=True)
        out_path = os.path.join(str(plot_dir),
                                f"{source_name}_to_{target_name}.png")

    if out_path is not None:
        fig.savefig(out_path, dpi=dpi, bbox_inches='tight')
        plt.close(fig)
        return out_path

    return fig

def plot_evaluate(source_fwhm, target_fwhm_v, D_v, Wm_v ):

    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(12,4))
    
    ax1.plot(target_fwhm_v, D_v, label='D', lw=4)
    ax1.set_xlabel("Gaussian FWHM")
    ax1.set_ylabel("D")
    #ax1.legend()
    ax2.plot(target_fwhm_v, Wm_v, label='W', lw=4)
    ax2.set_xlabel("Gaussian FWHM")
    ax2.set_ylabel(r"$W_{-}$")

    ax2.axhline(y=1, ls='--', c='k')
    ax2.axhline(y=0.5, ls='--', c='k')
    ax2.axhline(y=0.3, ls='--', c='k')

    out = np.interp(np.array([0.3, 0.5, 1.0]), Wm_v[::-1], target_fwhm_v[::-1])
    
    
    ax2.text(target_fwhm_v[-2]*0.8, 0.31, '{:.3f}"'.format(out[0]))
    ax2.text(target_fwhm_v[-2]*0.8, 0.51, '{:.3f}"'.format(out[1]))
    ax2.text(target_fwhm_v[-2]*0.8, 1.01, '{:.3f}"'.format(out[2]))
    
    out2 = np.interp(out,target_fwhm_v, D_v )
   
    for ii in range(len(out2)):
        ax1.axvline(x=out[ii], ls='--', c='k')
        ax1.text(out[ii], 0.11, '{:.3f}"'.format(out[ii]))


def find_safe_kernel(input_filter, detector_effects=True, save_kernels=True, verbose=False):

    # directories for PSF and kernels
    #psf_dir = '/'.join(path.dirname(path.realpath(jwst_kernels.__file__)).split('/')[:-2])+'/data/PSF/'
    # input_filter = {'camera':'MIRI', 'filter':'F2100W'}
    # input_filter = {'camera':'NIRCam', 'filter':'F300M'}
    # detector_effects =True
    # save_kernels=True
    # kernels_dir = '/'.join(path.dirname(path.realpath(jwst_kernels.__file__)).split('/')[:-2])+'/data/kernels/'
    source_psf, source_pixscale = read_PSF(input_filter, detector_effects=detector_effects)
    
    source_fwhm = fit_2d_gaussian(source_psf, pixscale=source_pixscale)
    #print('source FWHM', source_fwhm, source_pixscale)
    # Do a systematic search of the best Gaussian kernel by exploring kernels up to FWHM_Gauss = [1.05-2]*FWHM_source_PSF
    # Here we calcualted 11 kernels
    factor = np.linspace(1.05, 2, 18)
    target_fwhm_v = factor*source_fwhm
    size_kernel_asec = source_fwhm*10
    
    D_v, Wm_v = np.zeros(len(factor)), np.zeros(len(factor))
    #print(target_fwhm_v, size_kernel_asec)
    
    for ii, ff in enumerate(target_fwhm_v):
        if verbose==True:
            print('testing the nth PSF' +str(ii)+' with fwhm'+ '{:.3f}'.format(ff) )
        target_gaussian = {'fwhm': ff}
        kk = make_jwst_kernel_to_Gauss(input_filter, target_gaussian, 
                                       save_kernel=save_kernels, size_kernel_asec=size_kernel_asec, verbose=verbose)
        #plot_kernel(kk, save_plot=True, save_dir =None, want_convolve=True )
        D_v[ii], Wm_v[ii] = evaluate_kernel(kk)
        #print(D_v, Wm_v)
    
    Wfunct = interpolate.interp1d( Wm_v[::-1], target_fwhm_v[::-1])
    try:
        out = Wfunct(np.array([0.3, 0.5, 1.0]) )
    except:
        Warning('the output set of target PSFs attempted does not span the range in W=[0.3-1.0]. Try a larger range in the factor vector')
    
    print('Wm, very safe {:.3f}", safe {:.3f}", aggressive {:.3f}", source {:.3f}" '.format( 
        *out, source_fwhm))
    
    
    out2 = np.interp(out,target_fwhm_v, D_v )
    print('D, very safe {:.2f}, safe {:.2f}, aggressive {:.2f}'.format( *out2))
    
    outp = {'very_safe':out[0],'safe':out[1] , 'aggressive':out[2], 'source_fwhm':source_fwhm,
            'target_fwhm':target_fwhm_v , 'D_v':D_v, 'Wm_v':Wm_v }

    return(outp)

if __name__ == "__main__":
    input_filter = {'camera':'MIRI', 'filter':'F2100W'}
    out = find_safe_kernel(input_filter, save_kernels=True)
    
    plot_evaluate(out['source_fwhm'], out['target_fwhm'], out['D_v'], out['Wm_v'] )

