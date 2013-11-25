% test for lagrangian optimisation over a single square
%
%   Copyright (c) 2005 Gabriel Peyré

n = 32;

global wavelet_vm;
wavelet_vm = 0;

options.n = 512;    % size of the undelying big image
p = 16;
name = 'barb';
name = 'line';
disp('Loading image.');
options.use_wavelets = 1;
M = load_image_small_square(name, p, options);
y = l2error(M(:));

T = 10;

% factor 2 over-sampling
s = 6;
wavelet_vm = 0;
disp('Computing optimal direction.');
[MW1,theta] = compute_best_direction(M,T,s);
y1 = l2error(MW1(:));

% exhaustive search haar
s = Inf;
wavelet_vm = 0;
disp('Computing optimal direction.');
[MW2,theta] = compute_best_direction(M,T,s);
y2 = l2error(MW2(:));

% exhaustive search 7-9
s = Inf;
wavelet_vm = 4;
disp('Computing optimal direction.');
[MW3,theta] = compute_best_direction(M,T,s);
y3 = l2error(MW3(:));

q = 0.5*p^2;
sel = 1:q;
plot(   log2(sel), log2(y(sel)), ...
        log2(sel), log2(y1(sel)), ...
        log2(sel), log2(y2(sel)), ...
        log2(sel), log2(y3(sel)) );
legend('original', 'haarx2', 'haarxInf', '7/9xInf');
axis([0, log2(q), -2, max(log2(y))]);

% plot_geometry(theta, M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display of the Lagrangian
% compute tested directions
[Y,X] = meshgrid(0:n-1, 0:n-1); X = X(:); Y = Y(:);
X(1) = []; Y(1) = [];
Theta = atan2(Y(:),X(:));
Theta = unique(Theta);
Theta = [-Theta(end-1:-1:2); Theta]';
Theta = ( Theta + [Theta(2:end),Theta(1)+pi] )/2;
% compute Lagrangian
lambda = 3/(4*7);
% compute the lagrangian
L = []; % to store the 
for theta = Theta
    MW = perform_warped_wavelet(M,theta,1);
    % compute the quantized vector and number of bits
    [MWt,tmp,R] = perform_quantization(MW(:),T);
    % compute the approximation error
    E = sum( (MWt(:)-MW(:)).^2);
    % compute the lagrangian
    L = [L, E + lambda*R*t^2];
end

% the more the Lagrian as a frank minimum, the cleaner the geometry
% on the square is.
plot(Theta, L);
axis tight:
title('Lagrangian in function of geometry angle.');
ylabel('L(theta)');
xlabel('Theta')