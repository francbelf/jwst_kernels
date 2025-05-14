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
Install the code as a python package. Note that generating the PSFs on the fly requires the installation of stpfs (which is not installed automatically).

    python setup.py develop


Features
------------

Generates PSFs of the relevant bands using webbpsf seemlessly. Note that this package does **not** automatically install stpfs! You need to install it independently

Uses the Aniano+2011 algorithm to generate appropriate kernels for going between JWST bands and from a JWST band to a Gaussian.
Example usage to go between two JWST bands:

    from jwst_kernels.make_kernels import make_jwst_cross_kernel 

    input_filter = {'camera':'MIRI', 'filter':'F770W'}

    target_filter = {'camera':'MIRI', 'filter':'F2100W'}

    kk = make_jwst_cross_kernel(input_filter, target_filter)

Evaluate the kernels by finding the smallest safe Gaussian

    input_filter = {'camera':'NIRCam', 'filter':'F200W'}

    out = find_safe_kernel(input_filter, detector_effects=True) 

    print(out['safe'])

See examples for use of these functions in the example notebook <https://github.com/francbelf/jwst_kernels/blob/master/notebooks/examples.ipynb> 


.. _pyscaffold-notes:

Note
====

This project has been set up using PyScaffold 4.5. For details and usage
information on PyScaffold see https://pyscaffold.org/.
