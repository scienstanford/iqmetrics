% test for the quincunx wavelet transform

name = 'barb';
n = 512;
M = load_image(name, n);

Jmax = log2(n)-1;
Jmin = Jmax-2;

options.null = 0;
% Forward transform
fprintf('Computing forward transform ... ');
[MW,options.quincunx_filters] = perform_quincunx_wavelet_transform(M,Jmin,+1,options);
fprintf('done.\n');
% Backward transform
fprintf('Computing backward transform ... ');
M1 = perform_quincunx_wavelet_transform(MW,Jmin,-1,options);
fprintf('done.\n');

% error of reconstruction
disp([ 'Error of reconstruction ' num2str(psnr(M,M1), 4) 'dB.']);

% display
clf;
plot_quincunx_wavelet(MW, Jmin, options);

% select some sub-images
options.transform = 'quincunx';
m = 2*(Jmax-Jmin+1)+1;
nb = ceil(m/2);
k = 0;
clf;
for j=Jmax:-1:Jmin
    qq = 1:2;
    if j==Jmin
        qq(end+1) = 0;
    end
    for q=qq
        k = k+1;
        [selx,sely] = compute_quadrant_selection(j,q, options);
        subplot(2,nb,k);
        imagesc(MW(selx,sely));
        title(['j=' num2str(j) ', q=' num2str(q)]);
        axis image; axis off;
    end
end
colormap gray(256);