function plot_wavelet_quadtree(QT,Theta,j,q,M, plot_type)

% plot_wavelet_quadtree - plot a quadtree at a given wavelet level.
%
%   plot_wavelet_quadtree(QT,Theta,j,q,M, plot_type);
%
%   M is a wavelet transformed image.
%   QT and Theta are computed using the function 
%        [QT,Theta] = compute_wavelet_quadtree(M,Jmin,T,j_min,j_max,s);
%
%   Copyright (c) 2005 Gabriel Peyré

if nargin<6
    plot_type = 1;
end
if nargin<5
    M = [];
end

if q==1 % 1st quadrant
    selx = 1:2^j; sely = (2^j+1):2^(j+1);
elseif q==2
    selx = (2^j+1):2^(j+1); sely = 1:2^j;
else
    selx = (2^j+1):2^(j+1); sely = (2^j+1):2^(j+1);
end

plot_quadtree(QT(selx,sely),Theta(selx,sely), M(selx,sely), plot_type);