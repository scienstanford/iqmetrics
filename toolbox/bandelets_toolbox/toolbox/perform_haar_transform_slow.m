function y = perform_haar_transform_slow(x, Jmin, dir);

% perform_haar_transform - compute a wavelet haar transform
%
% y = perform_haar_transform(x, Jmin, dir);
%
%   This is the non-optimized Matlab version.
%   TODO : code inverse transform (the mex version does implement it).
%
%   Copyright (c) 2004 Gabriel Peyré

x = x(:);
n = length(x);
J = floor( log2(n) )-1;

if nargin<3
    dir=1;
end

if dir~=1 && dir~=-1
	error('dir should be either +1 or -1.');
end

if dir==-1
    warning('Inverse haar transform not implemented');    
    y = x;
    return;
end


% x contains the coarse scale signal
y = x;
for j=J:-1:0
    n = length(x);
    n1 = ceil(n/2);
    n2 = floor(n/2);
    
    % fine scale
    y(n1+1:n) = ( x(1:2:2*n2) - x(2:2:2*n2) )/sqrt(2);
    % coarse scale
    if n1==n2
        y(1:n1) = ( x(1:2:n) + x(2:2:n) )/sqrt(2);
    else
        y(1:n1-1) = ( x(1:2:n-1) + x(2:2:n) )/sqrt(2);
        y(n1) = x(n);
    end
    x = y(1:n1);
end

return;


% x contains the coarse scale signal
y = x;
for j=J:-1:0
    n = length(x);
    n1 = ceil(n/2);
    n2 = floor(n/2);
    
    if n1==n2
        % coarse scale
        y(1:n1) = ( x(1:2:n) + x(2:2:n) )/sqrt(2);
        % fine scale
        y(n1+1:n) = ( x(1:2:2*n2) - x(2:2:2*n2) )/sqrt(2);
        x = y(1:n1);
    else
        % coarse scale
        y(1:(n-3)/2) = ( x(1:2:n-3) + x(2:2:n-3) )/sqrt(2);
        y((n-1)/2) = ( x(n-2) + x(n-1) + x(n) )/sqrt(3);
        % fine scale
        y(n1:n-1) = ( x(1:2:n-1) - x(2:2:n-1) )/sqrt(2);
        y(n) = ( x(n-2) + x(n-1) - 2*x(n) )/sqrt(6);
        x = y(1:n1-1);
    end
    
end

return;