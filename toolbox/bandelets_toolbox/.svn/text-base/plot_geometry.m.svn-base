function plot_geometry(theta, M, str)

% plot_geometry - plot a geometry direction
%
%   plot_geometry(theta, M, str);
%
%   M should be a single dyadic square.
%   str is the style for the geometry (default str='r').
%
%   Copyright (c) 2005 Gabriel Peyré

if nargin<2
    M = [];
end
if nargin<3
    str = 'r';
end

hold on;

if ~isempty(M);
    n = size(M,1);
    x = linspace(0,1,n);
    imagesc(x,x,M);
end

pos = [0.5,0.5];
w = 1;
x = pos(1) + w/2*[cos(theta), -cos(theta)];
y = pos(2) + w/2*[sin(theta), -sin(theta)];
plot(x,y, str);

hold off;