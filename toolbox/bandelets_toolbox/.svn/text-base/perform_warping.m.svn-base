function [t,u,v,ts,us,vs,I] = perform_warping(M,theta,P)

% sampling location
[Y,X] = meshgrid(1:P,1:P);
% projection on orthogonal direction
t = -sin(theta)*X(:) + cos(theta)*Y(:);
u = cos(theta)*X(:) + sin(theta)*Y(:);
% order points in increasing order
[t,I] = sort(t);
u = u(I);
% sorted signal
v = M(I);
% remove small entries
J = find( abs(v)>1e-2 );
vs = v(J);
ts = t(J);
us = u(J);