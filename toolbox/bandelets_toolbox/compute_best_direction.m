function [MW,theta,L] = compute_best_direction(M,T,s)

% compute_best_direction - optimize the lagrangian over a single square
%
%   [MW,theta] = compute_best_direction(M,T,s);
%
%   M is a 2D image, s is a super-resolution factor
%       (the # of tested direction is 2*n*s)
%   T is the threshold (the higher, the most you want to compress)
%   theta is the optimal direction, ie. the one that gives minimizes the
%   Lagrangian
%       L(theta) = | M - F_theta(M) |^2 + m * T^2
%   where :
%       F_theta(M) = perform_warped_haar(M,theta,1)
%       m = #{ coefficients of |F_theta(M)| above T }
%
%   Copyright (c) 2005 Gabriel Peyré

if nargin<3
    s = 2;
end

n = size(M,1);

% samples the direction with a step of pi/(2*n*s)
if s~=Inf
    t = pi/(2*n*s);
    Theta = [t/2:t:pi-t/2, Inf];
else
    % perform exhaustive search
    [Y,X] = meshgrid(0:n-1, 0:n-1); X = X(:); Y = Y(:);
    X(1) = []; Y(1) = [];
    Theta = atan2(Y(:),X(:));
    Theta = unique(Theta);
    Theta = [-Theta(end-1:-1:2); Theta]';
    % take mid points
    Theta = ( Theta + [Theta(2:end),Theta(1)+pi] )/2;
    Theta = [Theta, Inf];
end

% Number of bit for coding geometry / no geometry
% You can try other coding strategy...
nbr_bits_geom = 1;
nbr_bits_nogeom = 1;

% conversion multiplier nbr_coefs<->nbr_bits
gamma0 = 7;
% lagrange multiplier
lambda = 3/(4*gamma0);
% number of bits for geometry
Rg = ceil( log2(length(Theta)) );

% compute the lagrangian
L = []; % to store the 
for theta = Theta
    MW = perform_warped_wavelet(M,theta,1);
    LT = [];
    for t = T(:)'
        % compute the quantized vector and number of bits
        [MWt,tmp] = perform_quantization(MW(:),t);
        % compute the number of coefficient above threshold
        nbr_coefs = sum( MW(:)>t );
        % estimate the number of bits needed to code the coefficients
        R = nbr_coefs*gamma0;
        % compute the approximation error 
        E = sum( (MWt(:)-MW(:)).^2);
        % add geometry bits
        if theta~=Inf
            R = R + nbr_bits_geom + Rg;
        else
            R = R + nbr_bits_nogeom;  % less bit for no geometry
        end
        % compute the lagrangian
        l = E + lambda*R*t^2;
        LT = [LT, l];
    end
    L = [L;LT];
end

% find minimum of lagrangian
[L,I] = min(L); L = L(:);
theta = Theta(I); theta = theta(:);
MW = perform_warped_wavelet(M,theta(1),1);