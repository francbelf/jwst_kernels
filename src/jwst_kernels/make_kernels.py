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
from jwst_kernels.kernel_core import MakeConvolutionKernel, profile, fit_2d_gaussian
from jwst_kernels.make_psf import makeGaussian_2D, read_PSF
from astropy.convolution import convolve


def make_jwst_cross_kernel(input_filter, target_filter, psf_dir=None, outdir=None,save_kernel=True,
                           common_pixscale=None, detector_effects=True,
                           naming_convention='PHANGS', verbose=False):
    '''Generates and saves the kernel necessary to convolve the image taken in a 
    JWST input filter into a JWST output filter. It works for both MIRI and NIRCam.
    

    Parameters
    ----------
    input_filter : dict
        Dictionary containing 'camera' and 'filter' keys
    target_filter : dict
        Dictionary containing 'camera' and 'filter' keys.
    psf_dir : str, optional
        Path to the directory where the JWST PSFs are saved. The default is None.
    outdir : str, optional
        Path to the directory where the kernels will be saved. The default is None.
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
    
    source_psf, source_pixscale = read_PSF(input_filter, detector_effects=detector_effects, psf_dir=psf_dir)

    target_psf, target_pixscale = read_PSF(target_filter, detector_effects=detector_effects, psf_dir=psf_dir)
    
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
                               verbose=verbose
                               )
    kk.make_convolution_kernel()
    dict_extension = {'DETEF': detector_effects}
    if outdir is None:
        outdir = '/'.join(path.dirname(path.realpath(jwst_kernels.__file__)).split('/')[:-2])+'/data/kernels/'
    if save_kernel==True:
        kk.write_out_kernel(outdir =outdir, add_keys =dict_extension ,naming_convention=naming_convention)
    return kk


     
def make_jwst_kernel_to_Gauss(input_filter, target_gaussian, psf_dir=None, outdir=None, save_kernel=True,
                          detector_effects=True,naming_convention='PHANGS', verbose=False, size_kernel_asec = None):
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

    source_psf, source_pixscale = read_PSF(input_filter, detector_effects=detector_effects, psf_dir=psf_dir)
    
    common_pixscale=source_pixscale
    target_pixscale= source_pixscale
    
    if size_kernel_asec is None:
        sz=int(target_gaussian['fwhm']/2.355*20/target_pixscale)
    else:
        sz=int(size_kernel_asec/target_pixscale)
    
    if sz%2==0:
        sz=sz+1

    #print(sz)
  
    yy, xx = np.meshgrid(np.arange(sz)-(sz-1)/2,np.arange(sz)-(sz-1)/2 )
  
    target_psf = makeGaussian_2D((xx, yy), (0,0), (
         target_gaussian['fwhm']/2.355/target_pixscale, \
             target_gaussian['fwhm']/2.355/target_pixscale) )
    target_name = 'gauss{:.2f}'.format(target_gaussian['fwhm'])
    #print(fit_2d_gaussian(target_psf, pixscale=source_pixscale))

    
    if verbose==True:
        print('making kernel from '+input_filter['filter']+' with source PSF to '+'Gaussan with FWHM {:.3f}'.format(target_gaussian['fwhm']))
    
    grid_size_arcsec = np.array([sz*target_pixscale, sz*target_pixscale])
   # print('grid_size_arcsec', grid_size_arcsec, grid_size_arcsec/target_pixscale)
    kk = MakeConvolutionKernel(source_psf=source_psf,
                               source_pixscale=source_pixscale,
                               source_name=input_filter['filter'],
                               target_psf=target_psf,
                               target_pixscale=target_pixscale,
                               target_name= target_name,
                               common_pixscale = common_pixscale,
                               grid_size_arcsec =grid_size_arcsec,
                               verbose=verbose)
    kk.make_convolution_kernel()
    dict_extension = {'DETEF': detector_effects}

    if outdir is None:
        outdir = '/'.join(path.dirname(path.realpath(jwst_kernels.__file__)).split('/')[:-2])+'/data/kernels/'
    if save_kernel==True:
        kk.write_out_kernel(outdir =outdir , add_keys =dict_extension, naming_convention=naming_convention)
    return kk
    
def plot_kernel(kk, save_plot=False, save_dir =None, want_convolve=True ):
    """Plots source and target PSF and the kernel
    

    Parameters
    ----------
    save_plot : default False
    save_dir: default None
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
        if save_dir is None:
            save_dir = '/'.join(path.dirname(path.realpath(jwst_kernels.__file__)).split('/')[:-2])+'/data/kernels/'
        
        plt.savefig(save_dir+kk.source_name+'_'+kk.target_name+'.png',
                    dpi=300)

if __name__ == "__main__":
    # example script
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
            
            kk = make_jwst_kernel_to_Gauss(input_filter, target_gaussian)
            print(kk.kernel.shape)
            plot_kernel(kk,save_plot=True, want_convolve=choic[k])

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
        
        kk = make_jwst_cross_kernel(input_filter, target_filter)
        plot_kernel(kk,save_plot=True)