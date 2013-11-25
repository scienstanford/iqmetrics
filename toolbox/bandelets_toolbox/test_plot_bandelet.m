% plot a bandelet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 16;
Jw = log2(64);     % wavelet scale  (2^Jw = size of the wavelet image)
Jb = 3;             % bandelet scale (1=finest)
k = 2;              % alpert function type
q = 3;              % wavelet direction
alpert_vm = [2,2];
theta = 1/sqrt(2);  % geometry orientation

str = ['bandelet_' 'Jb' num2str(Jb) '_Jw' num2str(Jw) '_k' num2str(k) '_q' num2str(q) '_n' num2str(n)];

options.use_mex = 0;



% scale localisation

% build a dirac function
[Y,X] = meshgrid(1:n,1:n);
t = -sin(theta)*X(:) + cos(theta)*Y(:);
pos = [t,X(:)]';

% test fwd just to get the info
[w,info] = perform_alpert_transform_2d(t,pos,alpert_vm, 1, options);
% select coefficient
I = find( info.k==k & info.l == max(info.l-Jb+1) );
u = I(round(end/2));

y = zeros(n^2,1);
y(u) = 1;

[w,info] = perform_alpert_transform_2d(y,pos,alpert_vm, -1, options);


M = reshape(w,n,n);
M = M+1e-5*(rand(n)-0.5);

I = find(M>=0);
M(I) = rescale(M(I));
I = find(M<0);
M(I) = -rescale(-M(I));


% interpolation
m = 1024;

[selx,sely] = compute_quadrant_selection(Jw,q);
selx = selx(end/2-n/2+1:end/2+n/2);
sely = sely(end/2-n/2+1:end/2+n/2);
MW = zeros(m,m);
MW(selx,sely) = M;

% backward transform
options.wavelet_vm = 3;
M1 = perform_wavelet_transform(MW,1,-1, options);

% extract approximatrly the support
w = 0.1;
a = 2^Jw/(2*m)+w;
A2 = M1( round(end/2-a*m+1):round(end/2+a*m), round(end/2-a*m+1):round(end/2+a*m) );

n1 = round( (1+w)*n/2 )*2;
n1 = n;
A1 = zeros(n1,n1);
A1( end/2-n/2+1:end/2+n/2, end/2-n/2+1:end/2+n/2 ) = M;

clf;
subplot(1,2,1);
imagesc(A1);
axis square; axis off;
colormap gray(256);

subplot(1,2,2);
imagesc(A2);
axis square; axis off;

% save
rep = 'images/';
warning off;
imwrite(rescale(A1), [rep str '_discrete.png'], 'png');
imwrite(rescale(A2), [rep str '_continuous.png'], 'png');
warning on;

return;



%%%
% plot a bandelet

global wavelet_vm;
wavelet_vm = 1;
theta = 1/sqrt(2);
n = 32;


% scale localisation
Jb = 1;

% build a dirac function
[Y,X] = meshgrid(1:n,1:n);
t = -sin(theta)*X(:) + cos(theta)*Y(:);
[tmp,I] = sort(t);


k = 3/2 * n/2^Jb;
y = zeros(n^2,1);
y(k) = 1;

if wavelet_vm==1
     y = reorder_coefs(y,1,-1);
end

Mb = zeros(n,n);
Mb(I) = y;

M = perform_warped_wavelet(Mb,theta,-1);
M = M+1e-5*(rand(n)-0.5);

I = find(M>=0);
M(I) = rescale(M(I));
I = find(M<0);
M(I) = -rescale(-M(I));

imagesc(M);
colormap gray(256);