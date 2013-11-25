toolbox_wavelets - wavelets related functions.

Orthogonal wavelet transform:
- perform_wavelet_transform - WAVELAB implementation of the wavelet transform.
Wavelet transform using lifting implementation (symmetric boundary condition):
- perform_lifting_transform / perform_lifting_transform_slow - 1D wavelet transform via lifting (general interface).
- perform_lifting_transform_byname - string based interface.
- perform_79_transform - biorthogonal 7/9 1D wavelet transform
- perform_wavelet_transform_isotropic - multidimensional isotropic (i.e. classical) wavelet transform.
- perform_wavelet_transform_hyperbolic - multidimensional hyperbolic (i.e. fully tensorial) wavelet transform.
Pyramid transform:
- perform_pyramid_transform - Laplacian-like pyramidal transform.
- perform_pyramid_transform_do - Minh Do Pyramidal transform (much better).
- perform_pyramid_transform_simoncelli - Steerable pyramid implementation of the Laplacian.
- perform_pyramid_transform_ti - translation-invariance pyramid (difference of Gaussian filterings).
Other transforms:
- perform_haar_transform - a simple but very fast 1D haar transform.
- perform_atrou_transform - compute the "a trou" wavelet transform, i.e. without subsampling (try to use either RWT or CWP2 when available).
- perform_cpx_dualtree_transform - complex dual tree transform.
- perform_steerable_transform - Steerable pyramid transform.
Compression and coding function:
- perform_jp2k_degradation - Perform JPEG2000 coding and decoding of wavelet coefficients.
- perform_jbig_coding.m - Perform JBig binary image coding (not related to wavelet).
- perform_spiht_coding - Perform SPIHT coding (slow) of the wavelet coefficients.
- perform_wavelet_arithmetic_coding - Perform simple arithmetic coding of the wavelet coefficients.
- perform_shannon_estimation_wavelets - compute the entropy of a wavelet transform.
Helpers functions:
- compute_quadrant_selection - compute the indices for selecting coefficients at a given scale and orientation.
- reorder_coefs - switch from inplace (results of lifting transform) to classical ordering.
- plot_wavelet - plot wavelet using Mallat's ordering.
- convert_wavelets2list - extract each sub-image.
Test scripts: see test_???.m files.

Installation note : you need to add the content of toolbox/ in your Matlab path.

Copyright (c) 2006 Gabriel Peyré