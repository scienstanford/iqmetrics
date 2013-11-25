iqmetrics
=========

Image Quality Metrics in Matlab

This is a repository for image quality metrics that have been published on the web.  Some of the scripts rely on functions that can be downloaded from isetBio (also a git repository). 

Most of the functions were compiled by Wen Lu in 2012 as a collection of "fidelity" or "reference" metrics. Reference metrics compare an original (ideal) and distorted (compressed, added noise, blurred, etc.) version of the original.  With the exception of SCIELAB, most metrics operate on RGB images.  We need to figure out whether these are assumed to be linear (gamma of 1.0) or sRGB (gamma of 2.2) images. 

We plan to add CPIQ metrics and compare predictions for various metrics. 


