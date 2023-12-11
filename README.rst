.. These are examples of badges you might want to add to your README:
   please update the URLs accordingly

    .. image:: https://api.cirrus-ci.com/github/<USER>/jwst_kernels.svg?branch=main
        :alt: Built Status
        :target: https://cirrus-ci.com/github/<USER>/jwst_kernels
    .. image:: https://readthedocs.org/projects/jwst_kernels/badge/?version=latest
        :alt: ReadTheDocs
        :target: https://jwst_kernels.readthedocs.io/en/stable/
    .. image:: https://img.shields.io/coveralls/github/<USER>/jwst_kernels/main.svg
        :alt: Coveralls
        :target: https://coveralls.io/r/<USER>/jwst_kernels
    .. image:: https://img.shields.io/pypi/v/jwst_kernels.svg
        :alt: PyPI-Server
        :target: https://pypi.org/project/jwst_kernels/
    .. image:: https://img.shields.io/conda/vn/conda-forge/jwst_kernels.svg
        :alt: Conda-Forge
        :target: https://anaconda.org/conda-forge/jwst_kernels
    .. image:: https://pepy.tech/badge/jwst_kernels/month
        :alt: Monthly Downloads
        :target: https://pepy.tech/project/jwst_kernels
    .. image:: https://img.shields.io/twitter/url/http/shields.io.svg?style=social&label=Twitter
        :alt: Twitter
        :target: https://twitter.com/jwst_kernels

.. image:: https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold
    :alt: Project generated with PyScaffold
    :target: https://pyscaffold.org/

|

============
jwst_kernels
============


   Make kernels to convolve JWST images taken in one band to the PSF of another, or to a Gaussian PSF of arbitrary FWHM.
   Implements the Aniano algorithm.


1.	Generate PSFs of the relevant bands
Install webbpsf and its associated data. 
IMPORTANT NOTE: installation of this package does not automatically install webbpsf! You need to install it independently
NOTE: latest version of webbpsf (v>1.0.0) are only available for python > 3.9. 
>> make_psf.py

2.	Generate the kernels
>> kernel_core.py 
fully general Aniano kernel generation script (from Tom with some modifications from me)
>>make_kernels.py 
uses classes implemented in kernel_core and applied them to two common cases:
	- jwst PSF to Gaussian (save_kernels_to_Gauss)
	- one JWST PSF to another (save_jwst_cross_kernel)
It only includes a function to plot the kernel and check the results after convolution with the original PSF. It now implements an option to consider detector effects.

3.	Assess the kernels
>> determine_safe_kernels.py
Looks in detail at one band to determine the “safe” Gaussian kernel according to Aninano



.. _pyscaffold-notes:

Note
====

This project has been set up using PyScaffold 4.5. For details and usage
information on PyScaffold see https://pyscaffold.org/.
