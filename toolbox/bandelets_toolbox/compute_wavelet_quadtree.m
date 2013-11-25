function [QT,Theta] = compute_wavelet_quadtree(M,Jmin,T,j_min,j_max,s, options)

% compute_wavelet_quadtree - compute the wavelet-bandelet quadtree
%
%   [QT,Theta] = compute_wavelet_quadtree(M,Jmin,T,j_min,j_max,s);
%
%   Copyright (c) 2005 Gabriel Peyr?

options.null = 1;
if nargin<2
    Jmin = 4;
end
if nargin<6
    s = Inf;
end
if nargin<5
    j_max = min(5, log2(size(M,1)));
end
if nargin<4
    j_min = 2;
end

if isfield(options, 'use_single_qt')
    use_single_qt = options.use_single_qt;
else
    use_single_qt = 0;
end

nT = length(T);

% perform the wavelet transform
MW = perform_wavelet_transform(M,Jmin, 1);

n = size(M,1); Jmax = 4;
QT = zeros(n,n,nT);
Theta = zeros(n,n,nT);

% compute the transform for each scale and each direction
for j=Jmax:-1:Jmin  % for each scale
    j_max = min(j_max, j);
    if use_single_qt
        MWc = zeros(2^j, 2^j, 3);
    end
    for q=1:3   % for each orientation
        [selx,sely] = compute_quadrant_selection(j,q);
%         disp(['--> computing quadtree at scale ' num2str(j) ' orientation ' num2str(q) '.']);
        if ~use_single_qt
            [QT(selx,sely,:),Theta(selx,sely,:)] = compute_quadtree(MW(selx,sely),T,j_min,j_max,s);
        else
            MWc(:,:,q) = MW(selx,sely);
        end
    end
    if use_single_qt
        [QTc,Thetac] = compute_quadtree(MWc,T,j_min,j_max,s);
        for q=1:3   % for each orientation
            [selx,sely] = compute_quadrant_selection(j,q);
            QT(selx,sely,:) = QTc;
            Theta(selx,sely,:) = Thetac;
        end
    end
end