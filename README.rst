.. These are examples of badges you might want to add to your README:
   please update the URLs accordingly
.. image:: https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold
    :alt: Project generated with PyScaffold
    :target: https://pyscaffold.org/

|

============
jwst_kernels
============


Make kernels to convolve JWST images taken in one band to the PSF of another, or to a Gaussian PSF of arbitrary FWHM.
Implements the Aniano algorithm.

Installation
------------

    python setup.py develop


Features
------------

Generates PSFs of the relevant bands using webbpsf seemlessly. Note that this package does **not** automatically install webbpsf! You need to install it independently
NOTE: latest version of webbpsf (v>1.0.0) are only available for python > 3.9. 

Uses the Aniano+2011 algorithm to generate appropriate kernels for going between JWST bands and from a JWST band to a Gaussian.
Example usage to go between two JWST bands:

    input_filter = {'camera':'MIRI', 'filter':'F770W'}

    target_filter = {'camera':'MIRI', 'filter':'F2100W'}

    kk = make_jwst_cross_kernel(input_filter, target_filter

Evaluate the kernels by finding the smallest safe Gaussian

    input_filter = {'camera':'NIRCam', 'filter':'F200W'}

    out = find_safe_kernel(input_filter, detector_effects=True) 
    
    print(out['safe'])

.. _pyscaffold-notes:

Note
====

This project has been set up using PyScaffold 4.5. For details and usage
information on PyScaffold see https://pyscaffold.org/.
