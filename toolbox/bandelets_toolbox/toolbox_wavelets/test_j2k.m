% test for jpeg2k compression
%
%   Copyright (c) 2005 Gabriel Peyré
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 512;
name = 'barb';
M = load_image(name, n);
Jmin = 3;

disp('Performing wavelet transform.');
MW = perform_wavelet_transform(M,Jmin,1);

bit_depth = 8;
true_comp_ratio = 12;

bpp = 0.5;
nbr_bits = (n^2)*bpp;

disp('Performing JPEG2k degradation.');
MW1 = perform_jp2k_degradation(MW,Jmin,nbr_bits,M);

% same but in 2 time
disp('Performing JPEG2k coding/decoding.');
options.Jmin = Jmin;
options.nbr_bits = nbr_bits;
% for i=1:20
stream = perform_wavelet_jp2k_coding(MW, options);
% end
MW2 = perform_wavelet_jp2k_coding(stream, options);
% should be zero
disp( ['Should be 0 : ' num2str(norme(MW1-MW2),2)] );

disp('Performing inverse wavelet transform.');
M1 = perform_wavelet_transform(MW1,Jmin, -1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display
ax(1) = subplot(1,2,1);
imagesc(M);
axis image; axis off;
title('Original');

ax(2) = subplot(1,2,2);
imagesc(M1);
axis image; axis off;
colormap gray(256);
title(['Compressed @' num2str(bpp, 2) 'bpp']);
linkaxes(ax,'xy');