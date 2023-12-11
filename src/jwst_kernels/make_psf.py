#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 19:54:42 2023

Use webbpsf to save a set of JWST PSF files for future use

@author: belfiore
"""

import numpy as np
from os import path
import jwst_kernels

import webbpsf


def makeGaussian_2D(X, M, S, normalise=False):
    gauss = np.exp(-np.power((X[0] - M[0])/S[0], 2.)/2)*np.exp(-np.power((X[1] - M[1])/S[1], 2.)/2)
    if normalise==True: gauss =gauss *1./(2.*np.pi*S[0]*S[1])
    return gauss


def save_miri_PSF(miri_psfs, output_dir='', **kwargs):
    """Generates MIRI PSF using webbpsf and saves then in output dir
    
    
    Parameters
    ----------
    miri_psfs : list
        list of MIRI filters.
        
    output_dir: string
        path to the output directory to save the PSF files
 
    """

    oversample_factor = kwargs.pop("oversample_factor", 4)
    detector_oversample = kwargs.pop("oversample_factor", 4)
    fov_arcsec = kwargs.pop("fov_arcsec", 19.98)
  
    for filter1 in miri_psfs:
        print('building PSF '+filter1)
    
        miri = webbpsf.MIRI()
        
        miri.filter = filter1
        psf_array = miri.calc_psf(oversample=oversample_factor,
               detector_oversample=detector_oversample,
               fov_arcsec=fov_arcsec, **kwargs)
        
        psf_array.writeto(output_dir+'MIRI_PSF_filter_'+miri.filter+'.fits',
                          overwrite=True)
    
def save_nircam_PSF(nircam_psfs, output_dir='', **kwargs):
    """Generates NIRCam PSF using webbpsf and saves then in output dir
    
    
    Parameters
    ----------
    nircam_psfs : list
        list of NIRCam filters.
        
    output_dir: string
        path to the output directory to save the PSF files
 
    """
    oversample_factor = kwargs.pop("oversample_factor", 4)
    detector_oversample = kwargs.pop("oversample_factor", 4)
    fov_arcsec = kwargs.pop("fov_arcsec", 10)
    
    for filter1 in nircam_psfs:
        print('building PSF '+filter1)
    
        nircam = webbpsf.NIRCam()
        
        nircam.filter = filter1
        psf_array = nircam.calc_psf(oversample=oversample_factor,
               detector_oversample=detector_oversample,
               fov_arcsec=fov_arcsec, **kwargs)
        
        psf_array.writeto(output_dir+'NIRCam_PSF_filter_'+nircam.filter+'.fits', 
                          overwrite=True)


if __name__ == "__main__":
    # example script
    # output directory where you want the JWST PSFs to be saved
    output_dir = '/'.join(path.dirname(path.realpath(jwst_kernels.__file__)).split('/')[:-2])+'/data/PSF/'

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

    for i in nircam_psfs:
        save_nircam_PSF(nircam_psfs, output_dir=output_dir)

    for i in miri_psfs:
         save_miri_PSF(miri_psfs, output_dir=output_dir)

