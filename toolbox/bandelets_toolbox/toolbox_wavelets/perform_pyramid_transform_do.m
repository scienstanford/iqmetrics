function y = perform_pyramid_transform_do(x, pfilt, k)

% perform_pyramid_transform_do - call the pyramidal transform of Do & Vetterli
%
%   y = perform_pyramid_transform_do(x, pfilt, k);
%
%   pfilt is the type of filter, can be '9-7', '5-3', 'Burt'
%
%   Use the algorihtm of 
%       Laplacian Pyramid Toolbox (version 1.0, June 2004)
%       http://www.ifp.uiuc.edu/~minhdo/software/
%   For more information see 
%       M. N. Do and M. Vetterli, Framing pyramids, 
%       IEEE Transactions on Signal Processing, 
%       vol. 51, pp. 2329-2342, Sep. 2003
%
%   Copyright (c) 2005 Gabriel Peyré

if nargin<2
    pfilt = '9-7';
end
if nargin<3
    k = 3;
end

if ~iscell(x)
    y = lpd(x, pfilt, k);
else
    y = lpr(x, pfilt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = lpd(x, pfilt, nlev)
% LPD   Multi-level Laplacian pyramid decomposition
%
%	y = lpdecn(x, pfilt, nlev)
%
% Input:
%   x:      input signal (of any dimension)
%   pfilt:  pyramid filter name (see PFILTERS)
%   nlev:   number of decomposition level
%
% Output:
%   y:      output in a cell vector from coarse to fine layers
%
% See also: LPR

% Get the pyramidal filters from the filter name
[h, g] = pfilters(pfilt);
% Decide extension mode
switch pfilt
    case {'9-7', '9/7', '5-3', '5/3', 'Burt'}
        extmod = 'sym';
    otherwise
        extmod = 'per';  
end
y = cell(1, nlev+1);
for n = 1:nlev
    [x, y{nlev-n+2}] = lpdec1(x, h, g, extmod);
end
y{1} = x;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = lpr(y, pfilt)
% LPR   Multi-level Laplacian pyramid reconstruction
%
%	x = lpr(y, pfilt)
%
% Input:
%   y:      output of a Laplacian pyramid in a cell vector
%   pfilt:  pyramid filter name (see PFILTERS)
%
% Output:
%   x:      reconstructed signal
%
% See also: LPD

% Get the pyramidal filters from the filter name
[h, g] = pfilters(pfilt);
% Decide extension mode
switch pfilt
    case {'9-7', '9/7', '5-3', '5/3', 'Burt'}
        extmod = 'sym';
    otherwise
        extmod = 'per';
end
x = y{1};
for n = 2:length(y)
    x = lprec1(x, y{n}, h, g, extmod);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c, d] = lpdec1(x, h, g, extmod)
% LPDEC1   One-level Laplacian pyramid decomposition
%
%	[c, d] = lpdec1(x, h, g)
%
% Input:
%   x:      input signal
%   h, g:   two biorthogonal 1-D lowpass filters
%   extmod: [optional] extension mode (default is 'per')
%
% Output:
%   c:      coarse signal at half size
%   d:      detail signal at full size
%
% See also:	LPREC1

if ~exist('extmod', 'var')
    extmod = 'per';
end
nd = ndims(x);
% Computer the coarse signal by filter and downsample
c = x;
for dim = 1:nd
    c = filtdn(c, h, dim, extmod, 0);
end
% Compute the detail signal by upsample, filter, and subtract
% Even size filter needs to be adjusted to obtain perfect reconstruction
adjust = mod(length(g) + 1, 2);
p = c;
for dim = 1:nd
    p = upfilt(p, g, dim, extmod, adjust);
end
d = x - p;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = lprec1(c, d, h, g, extmod)
% LPREC1   One-level Laplacian pyramid reconstruction
%
%	x = lprec1(c, d, h, g)
%
% Input:
%   c:      coarse signal at half size
%   d:      detail signal at full size
%   h, g:   two biorthogonal 1-D lowpass filters
%   extmod: [optional] extension mode (default is 'per')
%
% Output:
%   x:      reconstructed signal
%
% Note:     This uses a new reconstruction method by Do and Vetterli,
%           "Framming pyramids", IEEE Trans. on Sig Proc., Sep. 2003.
%
% See also:	LPDEC1

if ~exist('extmod', 'var')
    extmod = 'per';
end
nd = ndims(c);
% First, filter and downsample the detail image
r = d;
for dim = 1:nd
    r = filtdn(r, h, dim, extmod, 0);
end
% Then subtract the result from the coarse signal
p = c - r;
% Even size filter needs to be adjusted to obtain perfect reconstruction
adjust = mod(length(g) + 1, 2); 
% Then upsample and filter
for dim = 1:nd
    p = upfilt(p, g, dim, extmod, adjust);
end
% Final combination
x = p + d;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = filtdn(x, f, dim, extmod, shift)
% FILTDN   Filter and downsample (by 2) along a dimension
%
%       y = filtdn(x, f, dim, extmod, shift)
%
% Input:
%   x:      input signal
%   f:      1-D filter
%   dim:    the processing dimension
%   extmod: extension mode (e.g. 'per' or 'sym')
%   shift:  specifies the window over which filtering occurs
%
% Output:
%   y:      filtered and dowsampled signal
%
% Note:
%   The origin of the filter f is assumed to be floor(size(f)/2) + 1.
%   Amount of shift should be no more than floor((size(f)-1)/2).

% Skip singleton dimension
if size(x, dim) == 1    
    y = x;
    return
end
% Cell array of indexes for each dimension
nd = ndims(x);
I = cell(1, nd);
for d = 1:nd
    I{d} = 1:size(x,d);
end
% Border extend
n = size(x, dim);
hlf = (length(f) - 1) / 2;
% Amount of extension at two ends
e1 = floor(hlf) + shift;
e2 = ceil(hlf) - shift;
switch extmod
    case 'per'
        I{dim} = [ly-e1+1:n , 1:n , 1:e2];
    case 'sym'
        I{dim} = [e1+1:-1:2 , 1:n , n-1:-1:e2];
    otherwise
        error('Invalid input for EXTMOD')
end
y = x(I{:});
% Filter, downsample, and return only the 'valid' part
y = filter(f, 1, y, [], dim);
I{dim} = (1:2:n) + length(f) - 1;
y = y(I{:});









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = upfilt(x, f, dim, extmod, shift)
% UPFILT   Upsample (by 2) and filter along a dimension
%
%       y = upfilt(x, f, dim, extmod, shift)
%
% Input:
%   x:      input signal
%   f:      1-D filter
%   dim:    the processing dimension
%   extmod: extension mode (e.g. 'per' or 'sym')
%   shift:  specifies the window over which filtering occurs
%
% Output:
%   y:      upsampled and filtered signal
%
% Note:
%   The origin of the filter f is assumed to be floor(size(f)/2) + 1.
%   Amount of shift should be no more than floor((size(f)-1)/2).

% Skip singleton dimension
if size(x, dim) == 1
    y = x;
    return
end
nd = ndims(x);
% Consider column vectors as 1-D signals 
if nd == 2 & size(x, 2) == 1
    nd = 1;
end
% Cell array of indexes for each dimension
I = cell(1, ndims(x));
for d = 1:ndims(x)
    I{d} = 1:size(x,d);
end
% Upsample (by 2)
sx = size(x);
sx(dim) = 2*sx(dim);
y = zeros(sx);
I{dim} = 1:2:sx(dim);
y(I{:}) = x;
% Border extend
n = size(y, dim);
hlf = (length(f) - 1) / 2;
% Amount of extension at two ends
e1 = floor(hlf) + shift;
e2 = ceil(hlf) - shift;
switch extmod
    case 'per'
        I{dim} = [ly-e1+1:n , 1:n , 1:e2];  
    case 'sym'
        I{dim} = [e1+1:-1:2 , 1:n , n-1:-1:e2];
    otherwise
        error('Invalid input for EXTMOD')
end
y = y(I{:});
% Filter and return only the 'valid' part
y = filter(f, 1, y, [], dim);
I{dim} = (1:n) + length(f) - 1;
y = y(I{:});






function [h, g] = pfilters(fname)
% PFILTERS    Generate filters for the Laplacian pyramid
%
%	[h, g] = pfilters(fname)
%
% Input:
%   fname:  Name of the filters, including the famous '9-7' filters
%           and all other available from WFILTERS in Wavelet toolbox
%
% Output:
%   h, g:   1D filters (lowpass for analysis and synthesis, respectively)
%           for seperable pyramid
switch fname
    case {'9-7', '9/7'}
        h = [.037828455506995 -.023849465019380 -.11062440441842 ...
             .37740285561265];	
        h = [h, .85269867900940, fliplr(h)];
    
        g = [-.064538882628938 -.040689417609558 .41809227322221];
        g = [g, .78848561640566, fliplr(g)];
    case {'5-3', '5/3'}
        h = [-1, 2, 6, 2, -1] / (4 * sqrt(2));
        g = [1, 2, 1] / (2 * sqrt(2));
        
    case {'Burt'}
        h = [0.6, 0.25, -0.05];
        h = sqrt(2) * [h(end:-1:2), h];
	
        g = [17/28, 73/280, -3/56, -3/280];
        g = sqrt(2) * [g(end:-1:2), g];
    otherwise
        [h, g] = wfilters(fname, 'l');
end