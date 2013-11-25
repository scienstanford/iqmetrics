function [MB,m_geom] = perform_wavelet_bandelet_transform(M,Jmin,QT,Theta,dir, options)

% perform_wavelet_bandelet_transform - perform the wavelet-bandelet transform
%
%   [MB,m_geom] = perform_wavelet_bandelet_transform(M,Jmin,QT,Theta,dir);
%
%   Copyright (c) 2005 Gabriel Peyr?

options.null = 1;
if nargin<2
    Jmin = 4;
end

if isfield(options, 'use_single_qt')
    use_single_qt = options.use_single_qt;
else
    use_single_qt = 0;
end

if size(QT,3)>1 || size(Theta,3)>1
    warning('Quadtree seems to have been computed using multiple thresholds.');
end

MB = M;
if dir==1
    % perform the wavelet transform
    MB = perform_wavelet_transform(M,Jmin, 1);
end

n = size(M,1); Jmax = 4; 
m_geom = 0;
% compute the transform for each scale and each direction
for j=Jmax:-1:Jmin  % for each scale
    for q=1:3   % for each orientation
        if q==1 % 1st quadrant
            selx = 1:2^j; sely = (2^j+1):2^(j+1);
        elseif q==2
            selx = (2^j+1):2^(j+1); sely = 1:2^j;
        else
            selx = (2^j+1):2^(j+1); sely = (2^j+1):2^(j+1);
        end
        [MB(selx,sely),m] = perform_bandelet_transform(MB(selx,sely),QT(selx,sely),Theta(selx,sely),dir);
        if use_single_qt
            m = m/3;
        end
        m_geom = m_geom + m;
    end    
end

if dir==-1
    % perform the inverse wavelet transform
    MB = perform_wavelet_transform(MB,Jmin, -1);
end