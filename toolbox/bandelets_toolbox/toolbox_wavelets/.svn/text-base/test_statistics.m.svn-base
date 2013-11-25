% some statistical tests


imrep = 'data/';
imrep = '';

name = 'hair';
name = 'line_vertical_blurred';
name = 'sawtooth';
name = 'bwcircle';
name = 'line_vertical';
name = 'flow2';
name = 'lena';
name = 'flow2';
name = 'barb';
name = 'noise';


rep = ['results_stats/' name '/'];

save_images = 1;

if exist(rep)~=7
    mkdir(rep);
end

n = 512;
if strcmp(name, 'noise')
    M = 30*randn(n);
else
    M = load_image([imrep name], n, options);
end


options.wavelet_type = 'biorthogonal';
options.wavelet_vm = 4;
Jmin = 4;
Jmax = log2(n)-1;

MW = perform_wavelet_transform(M, Jmin, +1, options);

j = Jmax; q = 1;
[selx,sely] = compute_quadrant_selection(j,q);
MWs = MW(selx,sely);
MWs = MW;

if save_images
    warning off;
    imwrite(rescale(M),[rep name 'original.png'], 'png');
    warning off;
end

% compute histogram
nb_bins = 501;
[N,X] = histo(MWs(:),nb_bins);
N = N/prod( size(MWs) );
plot( X, log2(N) );
axis( [-400 400 -16 max(log2(N))] );

if save_images
    saveas(gcf, [rep name '_wav_hiso.eps'], 'epsc');
    saveas(gcf, [rep name '_wav_hiso.png'], 'png');
end

x = [];
x_cond = [];
q = 1;
dist = 4;
for j=Jmax:-1:Jmin
    [selx,sely] = compute_quadrant_selection(j,q);
    MWs = MW(selx,sely);
    a = MWs(1:end-dist,:); 
    b = MWs(1+dist:end,:); 
    x = [x; a(:)];
    x_cond = [x_cond; b(:)];
end

p = 20;
pc = 20;
T = 3;
Tc = T*p/pc;
[tmp, xt] = perform_quantization(x, Tc, 1);
[tmp, xt_cond] = perform_quantization(x_cond, Tc, 1);

H = compute_symmetric_conditional_histogram(xt, xt_cond, p, pc);

J = find(H~=0);
I = find(H==0);
H(I) = min(H(J));

imagesc([-p*T, T*p], [-pc*Tc, Tc*pc], log2(H));
ylabel('Neighbor');
xlabel('Pixel');
colormap gray(256);
axis image; axis xy;

if save_images
    saveas(gcf, [rep name '_wav_hisocond.eps'], 'epsc');
    saveas(gcf, [rep name '_wav_hisocond.png'], 'png');
end


return;

% compute flow transform
options.nbr_scales = 4;
options.search_width = 2;
options.windows_width = 4;
options.enforce_direction = 1;
options.use_mex=1;
[Theta,MB,M1,F,ThetaAbs] = compute_bandelet_flow(M, options);