function distortion = WNISM_RR(origImg,distImg)

% A simple test for the reduced-reference image quality
% assessment algorithm published in
% Z. Wang and E. P. Simoncelli, "Reduced-reference image
% quality assessment using a wavelet-domain natural image
% statistic model," Human Vision and Electronic Imaging X,
% Prof. SPIE, vol. 5666, San Jose, CA, Jan. 2005

% set path
path(path, 'RRIQA/sender_side');
path(path, 'RRIQA/receiver_side');
path(path, 'RRIQA/common_to_both_sides');
	% read reference image
features_in_bits = sender_feature_extraction(origImg);

% receiver side: feature extraction and distortion measurement
	% read distorted image
distortion = receiver_distortion_measure(distImg, features_in_bits);