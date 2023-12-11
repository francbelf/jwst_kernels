#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 14:07:52 2022

This script uses webbpsf to generate a fresh version of the JWST kernels
Useful in the future if webbpsf is updated

@author: belfiore
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import table
import copy
from os import path

import jwst_kernels
from jwst_kernels.kernel_core import MakeConvolutionKernel, profile
from jwst_kernels.make_psf import save_miri_PSF, save_nircam_PSF, makeGaussian_2D
from astropy.convolution import convolve

def make_jwst_cross_kernel(input_filter, target_filter, psf_dir='', outdir='',
                           common_pixscale=None, detector_effects=True,
                           naming_convention='PHANGS'):
    '''Generates and saves the kernel necessary to convolve the image taken in a 
    JWST input filter into a JWST output filter. It works for both MIRI and NIRCam.
    

    Parameters
    ----------
    input_filter : dict
        Dictionary containing 'camera' and 'filter' keys
    target_filter : dict
        Dictionary containing 'camera' and 'filter' keys.
    psf_dir : str, optional
        Path to the directory where the JWST PSFs are saved. The default is ''.
    outdir : str, optional
        Path to the directory where the kernels will be saved. The default is ''.
    detector_effects: bool, default: True
        Whether to include detector effects in the JWST PSFs when generating kernels

    Returns
    -------
    kk : MakeConvolutionKernel
        Object containing the kernel.
        
    Notes
    -------
    If the necessary JWST PSF is not found in psf_dir, the code will use webbpsf
    to generate the PSF. This requires webbpsf to be installed and the necessary
    files to have been added to the path. For more details see 
    https://webbpsf.readthedocs.io/en/latest/installation.html
    '''
    
    if detector_effects ==False:
        extension = "PRIMARY"
    else:
        extension = 'OVERDIST'
  

    source_psf_path = psf_dir+input_filter['camera']+'_PSF_filter_'+\
            input_filter['filter']+'.fits'
    try:
        source_psf = fits.open(source_psf_path)

    except FileNotFoundError:
        print('generating PSF with webbpsf!')
        if input_filter['camera']=='MIRI':
            save_miri_PSF([input_filter['filter']], output_dir=psf_dir)
            source_psf = fits.open(source_psf_path)
            
        if input_filter['camera']=='NIRCam':
            save_nircam_PSF([input_filter['filter']], output_dir=psf_dir)
            source_psf = fits.open(source_psf_path)
    
    source_psf=source_psf[extension]
    source_pixscale = source_psf.header['PIXELSCL']
            
    target_psf_path = psf_dir+target_filter['camera']+'_PSF_filter_'+\
            target_filter['filter']+'.fits'
    try:
        target_psf = fits.open(target_psf_path)
       
    except FileNotFoundError:
        print('generating PSF with webbpsf!')
        if input_filter['camera']=='MIRI':
            save_miri_PSF([target_filter['filter']], output_dir=psf_dir)
            target_psf = fits.open(target_psf_path)
            
        if input_filter['camera']=='NIRCam':
            save_nircam_PSF([target_filter['filter']], output_dir=psf_dir)
            target_psf = fits.open(target_psf_path)
    
    target_psf=target_psf[extension]
    target_pixscale = target_psf.header['PIXELSCL']
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
                               common_pixscale = common_pixscale,
                               grid_size_arcsec =grid_size_arcsec,
                               verbose=True
                               )
    kk.make_convolution_kernel()
    dict_extension = {'extension': extension}
    kk.write_out_kernel(outdir =outdir, add_keys =dict_extension ,naming_convention='PHANGS')
    return kk


     
def save_kernels_to_Gauss(input_filter, target_gaussian, psf_dir='', outdir='', 
                          detector_effects=False,naming_convention='PHANGS'):
    '''Generates and saves the kernel necessary to convolve the image taken in a 
    JWST input filter into a JWST output filter. It works for both MIRI and NIRCam.
    

    Parameters
    ----------
    input_filter : dict
        Dictionary containing 'camera' and 'filter' keys
    target_gaussian : dict
        Dictionary containing a 'fwhm' key.
    psf_dir : str, optional
        Path to the directory where the JWST PSFs are saved. The default is ''.
    outdir : str, optional
        Path to the directory where the kernels will be saved. The default is ''.
    detector_effects: bool, default: True
        Whether to include detector effects in the JWST PSFs when generating kernels

    Returns
    -------
    kk : MakeConvolutionKernel
        Object containing the kernel.
        
    Notes
    -------
    If the necessary JWST PSF is not found in psf_dir, the code will use webbpsf
    to generate the PSF. This requires webbpsf to be installed and the necessary
    files to have been added to the path. For more details see 
    https://webbpsf.readthedocs.io/en/latest/installation.html
    '''

    source_psf_path = psf_dir+input_filter['camera']+'_PSF_filter_'+\
            input_filter['filter']+'.fits'
    if detector_effects ==False:
        extension = "PRIMARY"
    else:
        extension = 'OVERDIST'
    
    source_psf_path = psf_dir+input_filter['camera']+'_PSF_filter_'+\
            input_filter['filter']+'.fits'
    try:
        source_psf = fits.open(source_psf_path)

    except FileNotFoundError:
        print('generating PSF with webbpsf!')
        if input_filter['camera']=='MIRI':
            save_miri_PSF([input_filter['filter']], output_dir=psf_dir)
            source_psf = fits.open(source_psf_path)
            
        if input_filter['camera']=='NIRCam':
            save_nircam_PSF([input_filter['filter']], output_dir=psf_dir)
            source_psf = fits.open(source_psf_path)
    
    source_psf=source_psf[extension]
    source_pixscale = source_psf.header['PIXELSCL']
    
    common_pixscale=source_pixscale
    target_pixscale= source_pixscale

    sz=int(target_gaussian['fwhm']/2.355*10/target_pixscale)
    if sz%2==0:
        sz=sz+1
  
    yy, xx = np.meshgrid(np.arange(sz)-(sz-1)/2,np.arange(sz)-(sz-1)/2 )
  
    target_psf = makeGaussian_2D((xx, yy), (0,0), (
         target_gaussian['fwhm']/2.355/target_pixscale, \
             target_gaussian['fwhm']/2.355/target_pixscale) )
    target_name = 'gauss{:.2f}'.format(target_gaussian['fwhm'])
    
    grid_size_arcsec = np.array([sz*target_pixscale, sz*target_pixscale])
    print('grid_size_arcsec', grid_size_arcsec/target_pixscale)
    kk = MakeConvolutionKernel(source_psf=source_psf,
                               source_pixscale=source_pixscale,
                               source_name=input_filter['filter'],
                               target_psf=target_psf,
                               target_pixscale=target_pixscale,
                               target_name= target_name,
                               common_pixscale = common_pixscale,
                               grid_size_arcsec =grid_size_arcsec,
                               verbose=True)
    kk.make_convolution_kernel()
    dict_extension = {'extension': extension}
    kk.write_out_kernel(outdir =outdir , add_keys =dict_extension, naming_convention='PHANGS')
    return kk
        
# def get_copt_fwhm(gal_name):
#     """For a given PHANGS galaxy, get the FWHM of the copt MUSE data
#     """

#     t= table.Table.read('muse_dr2_v1.fits')
#     ii = t['name']==gal_name
#     copt_fwhm = float(t[ii]['muse_copt_FWHM'])
#     return copt_fwhm
    
    
def plot_kernel(kk, save_plot=False, save_dir ='', want_convolve=True ):
    """Plots source and target PSF and the kernel
    

    Parameters
    ----------
    save_plot : default False
    save_dir: default ''
    want_convolve: default True
        Whethere the plot should show the result of convolving the kernel with the soruce PSF (to visually check the 
        goodness of the kernel). It can get slow for large kernels.

    Returns
    -------
    None.

    """
    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(12,4))

    ax1.imshow(np.log10(kk.source_psf/np.max(kk.source_psf)), vmax=0, vmin=-4);
    ax1.set_title(kk.source_name)
    
    ax2.imshow(np.log10(kk.target_psf/np.max(kk.target_psf)), vmax=0, vmin=-4);
    ax2.set_title(kk.target_name)
    
    extent = int(10*kk.target_fwhm/kk.common_pixscale/2)
    ax3.plot(*profile(kk.source_psf/np.max(kk.source_psf), 
                      bins=np.linspace(0, 6*kk.target_fwhm,extent) ,
                      pixscale=kk.common_pixscale), 
             c='b', label=kk.source_name);
    ax3.plot(*profile(kk.target_psf/np.max(kk.target_psf), 
                      bins=np.linspace(0, 6*kk.target_fwhm, extent),
                      pixscale=kk.common_pixscale),
             c='k', label=kk.target_name, lw=5);
    
    ax3.plot(*profile(kk.kernel/np.max(kk.kernel), 
                      bins=np.linspace(0, 6*kk.target_fwhm, extent),
                      pixscale=kk.common_pixscale),
             c='g', label='kernel');
    if want_convolve==True:
        target_conv = convolve(kk.source_psf, kk.kernel)
        ax3.plot(*profile(target_conv/np.max(target_conv), 
                        bins=np.linspace(0, 6*kk.target_fwhm, extent),
                        pixscale=kk.common_pixscale),
                c='r', label='model', ls='-')
        
        xx, targetprof = profile( kk.target_psf/np.max(kk.target_psf), bins=np.linspace(0, 6*kk.target_fwhm, extent),
                        pixscale=kk.common_pixscale)
        trash, convprof = profile( target_conv/np.max(target_conv), bins=np.linspace(0, 6*kk.target_fwhm, extent),
                        pixscale=kk.common_pixscale)
        ax3.plot(xx, (targetprof-convprof),
                c='r', label='residual', ls='--')
    ax3.legend()
    #ax3.set_yscale('log')
    ax3.set_ylim([-0.1, 1.1])
    ax3.set_xlim([0,6*kk.target_fwhm])
    if save_plot ==True:
        plt.savefig(save_dir+kk.source_name+'_'+kk.target_name+'.png',
                    dpi=300)

if __name__ == "__main__":
    # example script
    # output directory where you want the JWST PSFs to be saved
    output_dir = '/'.join(path.dirname(path.realpath(jwst_kernels.__file__)).split('/')[:-2])+'/data/PSF/'
    output_dir_kernels = '/'.join(path.dirname(path.realpath(jwst_kernels.__file__)).split('/')[:-2])+'/data/kernels/'

    #list of the PHANGS-JWST filters, others can be added if necessary
    nircam_psfs = [
        'F150W',
        'F187N', 
        'F200W',
        'F300M',
        'F335M',
        'F360M',
    ]

    miri_psfs = [
        'F770W',
        'F1000W',
        'F1130W',
        'F2100W',
    ]

    #Generate kernels going from all the JWST filters into Gaussian PSF of specific FWHM
    all_PSFs = nircam_psfs+ miri_psfs
    all_cameras = ['NIRCam']*len(nircam_psfs) + ['MIRI']*len(miri_psfs)
    gauss_fwhm = [0.85, 0.9, 1.0, 4., 7.5, 15] # list of target Gaussian FWHM
    choic = [True, True, True, False, False, False] # turn off the convolution when plotting for larger kernels

    for k, jj in enumerate(gauss_fwhm):
        for ii in range(len(all_PSFs)):
            print( all_cameras[ii], all_PSFs[ii], ' to Gauss ' + '{:.3f}'.format(jj))
            input_filter = {'camera':all_cameras[ii], 'filter':all_PSFs[ii]}
            target_gaussian = {'fwhm':jj}
            
            kk = save_kernels_to_Gauss(input_filter, target_gaussian,
                                        psf_dir=output_dir, outdir=output_dir_kernels,
                                        detector_effects=True, naming_convention='PHANGS')
            print(kk.kernel.shape)
            plot_kernel(kk,save_plot=True, save_dir=output_dir_kernels, want_convolve=choic[k])

        # %%
    #Generate kernels going from all the JWST filters into F2100W 
    # (except F2100W--> F2100W which would be silly of course)
    miri_psfs_n = copy.copy(miri_psfs)
    miri_psfs_n.remove('F2100W')
    all_PSFs = nircam_psfs+ miri_psfs_n
    all_cameras = ['NIRCam']*len(nircam_psfs) + ['MIRI']*len(miri_psfs_n)

    for ii in range(len(all_PSFs)):
        print( all_cameras[ii], all_PSFs[ii], ' to F2100W')
        input_filter = {'camera':all_cameras[ii], 'filter':all_PSFs[ii]}
        target_filter = {'camera':'MIRI', 'filter':'F2100W'}
        
        kk = save_jwst_cross_kernel(input_filter, target_filter,
                                    psf_dir=output_dir, outdir=output_dir_kernels, naming_convention='PHANGS')
        plot_kernel(kk,save_plot=True, save_dir=output_dir_kernels)