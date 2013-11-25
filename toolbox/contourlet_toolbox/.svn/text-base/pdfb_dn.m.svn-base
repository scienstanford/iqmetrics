function xdn = pdfb_dn(x, pfilt, dfilt, nlevs, nstd, th, sigma)
% PDFB_DN   Image denoising using the contourlet or PDFB transform
%
%   xdn = pdfb_dn(x, pfilt, dfilt, nlevs, nstd, [th, sigma])
%
% Input:
%   x:      noisy input image
%   pfilt:  filter name for the pyramidal decomposition step
%   dfilt:  filter name for the directional decomposition step
%   nlevs:  vector of numbers of directional filter bank decomposition levels 
%   nstd:   noise standard deviation in the PDFB domain.  This is computed by:
%           nstd = pdfb_nest(size(x, 1), size(x, 2), pfilt, dfilt, nlevs);
%   th:     [optional] scale for threshold; default 3 for 3*sigma thresholding
%   sigma:  [optional] noise standard deviation in the image domain
%
% Output:
%   xdn:    denoised image

if ~exist('th', 'var')
    th = 3;
end

% Contourlet transform
y = pdfbdec(x, pfilt, dfilt, nlevs);
[c, s] = pdfb2vec(y);

% Estimage standard deviation of additive Gaussian white noise if not given
if ~exist('sigma', 'var');
    return;    
end

cth = th * sigma * sqrt(nvar);

% Slightly different thresholds for the finest scale
fs = s(end, 1);
fssize = sum(prod(s(find(s(:, 1) == fs), 3:4), 2));
cth(end-fssize+1:end) = (4/3) * cth(end-fssize+1:end);

c = c .* (abs(c) > cth);

% Reconstruction
y = vec2pdfb(c, s);
xdn = pdfbrec(y, pfilt, dfilt);

