#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 18:10:43 2023

@author: belfiore
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits, ascii
import astropy.table as table
from jwst_kernels.kernel_core import  fit_2d_gaussian
from jwst_kernels.make_kernels import make_jwst_kernel_to_Gauss, plot_kernel
from jwst_kernels.make_psf import  read_PSF
from astropy.convolution import convolve
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
    # Do a systematic search of the best Gaussian kernel by exploring kernels up to FWHM_Gauss = [1.04-1.5]*FWHM_source_PSF
    # Here we calcualted 11 kernels
    factor = np.linspace(1.05, 2, 18)
    target_fwhm_v = factor*source_fwhm
    size_kernel_asec = source_fwhm*10
    
    D_v, Wm_v = np.zeros(len(factor)), np.zeros(len(factor))
    #print(target_fwhm_v, size_kernel_asec)
    
    for ii, ff in enumerate(target_fwhm_v):
        if verbose==True:
            print('testing the nth PSF' +str(ii)+' with fwhm'+ '{:.f3}'.format(ff) )
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

#%%
