function [QT,Theta] = compute_quadtree(M,T,j_min,j_max,s)

% bandelet_transform_fwd - compute the quadtree that optimize the Lagrangian.
%
%   [QT,Theta] = compute_quadtree(M,T,j_min,j_max,s);
%
%   M is the original image
%   T is the selected threshold (the higher, the most compressed the data)
%   j_min is the depth minimum of the QT 
%       (ie 2^j_min is the size minimum of the square).
%   j_max is the depth maximum of the QT        [default : min(5,log2(n))]
%       (ie 2^j_max is the size maximum of the square).
%   s is the super-resolution for the geometry [default 2]
%
%   QT is an image representing the levels of the quadtree.
%   Theta is an image representing the optimal angle choosed on each
%       quadtree (Inf token for no geometry).
%
%   Copyright (c) 2005 Gabriel Peyr?

if nargin<5
    s = Inf;
end
if nargin<4
    j_max = min(5, log2(size(M,1)));
end
if nargin<3
    j_min = 2;
end

n = size(M,1);

% number of thresholds
nT = length(T);

% Lagrangian multiplicator
lambda = 3/(4*7);

QT = zeros(n,n,nT)+j_min;
Theta = zeros(n,n,nT);
% L = zeros(n/2^j_min,n/2^j_min,nT);   % the current lagrangian
L = zeros(floor(n/2^j_min),floor(n/2^j_min),nT);   % bug: n/2^j_min and n/2^j_min maybe not integer

% first compute bandelet approximation for each square
for kx=0:n/2^j_min-1
    for ky=0:n/2^j_min-1
        selx = kx*2^j_min+1:(kx+1)*2^j_min;
        sely = ky*2^j_min+1:(ky+1)*2^j_min;
        % compute the optimal direction on this square
        [tmp,theta,l] = compute_best_direction(M(selx,sely,:),T,s);
        for i=1:nT
            Theta(selx,sely,i) = theta(i);
            L(kx+1,ky+1,i) = l(i);
        end
    end
end

% perform the bottom-up procedure, trying to merge 
% 4 small squares into 1 big square if it decreases the 
% lagrangian
for j=j_min+1:j_max
    L1 = zeros(n/2^j);  % new lagrangian for this size of squares
    for kx=0:n/2^j-1
        for ky=0:n/2^j-1
            selx = kx*2^j+1:(kx+1)*2^j;
            sely = ky*2^j+1:(ky+1)*2^j;
            % the lagrangian of the 4 squares splited (add the gamma 
            % penalty because of the additional split)
            l_sum = L(2*kx+1,2*ky+1,:) + L(2*kx+2,2*ky+1,:) + L(2*kx+1,2*ky+2,:) + L(2*kx+2,2*ky+2,:) + lambda*reshape(T,1,1,nT).^2;
            l_sum = l_sum(:);
            % the lagrangian of the 4 squares once merged
            [tmp,theta,l] = compute_best_direction(M(selx,sely,:),T,s);
            % find the indices where l<l_sum, ie we must merge
            I = find(l<l_sum);
            J = find(l>=l_sum);
            % MERGE
            for i=I'
              L1(kx+1,ky+1,i) = l(i);
              QT(selx,sely,i) = j;
              Theta(selx,sely,i) = theta(i);
            end
            % KEEP
            for i=J'
              L1(kx+1,ky+1,i) = l_sum(i);
            end
        end
    end
    L = L1;
end