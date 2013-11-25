function M = perform_warped_wavelet(M,theta,dir)

% perform_warped_wavelet - perform a warped haar transform
%
%   MW = perform_warped_wavelet(M,theta,dir);
%
%   Copyright (c) 2005 Gabriel Peyré

if theta==Inf   % special token : no geometry
    return; % nothing to do
end

if nargin<3
    dir = 1;
end

if size(M,3)>1
    for i=1:size(M,3)
        M(:,:,i) = perform_warped_wavelet(M(:,:,i),theta,dir);
    end
    return;
end

n = size(M,1);
% sampling location
[Y,X] = meshgrid(1:n,1:n);
% projection on orthogonal direction
t = -sin(theta)*X(:) + cos(theta)*Y(:);
% order points in increasing order
[tmp,I] = sort(t);
Jmin = 1;
M(I) = perform_wavelet_transform(M(I),Jmin,dir);