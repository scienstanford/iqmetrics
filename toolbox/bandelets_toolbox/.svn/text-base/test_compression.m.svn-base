% test the quadtree construction in wavelet domain
%
%   Copyright (c) 2005 Gabriel Peyré

% set to 0 if you don't want to save the images
save_images = 0;

name = 'polygons_blurred';
name = 'barb';
rep = 'images/';

if save_images
    warning off;
    mkdir(rep);
    warning on;
end

global wavelet_vm;
wavelet_vm = 0;

n = 512;
s = Inf;    % super-resolution on the geometry

% set to:
%   1 if you want the 3 wavelets bands H/V/D to share the same quadtree
%   0 for 3 different quadtree (more adaptivity but more bits for geometry).
options.use_single_qt = 0;

M = load_image(name, n);

% reduce size for speed up
n = 128;
M = M(end/2-n/2+1:end/2+n/2, end/2-n/2+1:end/2+n/2);

% remove high transition by projecting
% on coarse scale wavelet domain
MW = perform_wavelet_transform(M,log2(n),1);
M = MW(1:n,1:n);

if save_images
    warning off;
    imwrite( M/256, [rep 'original.png'], 'png' )
    warning on;
end

Jmin = 4;
Jmax = log2(n)-1;
j_min = 2;
j_max = 4;
T = [1:0.5:5, 7:2:15, 20];
disp('Computing quadtree.');
if ~exist('QT')
    [QT,Theta] = compute_wavelet_quadtree(M,Jmin,T,j_min,j_max,s, options);
end

MW = perform_wavelet_transform(M,Jmin, 1);

% to reccord the results
psnr_wav = [];
bit_wav = [];
psnr_band = [];
bit_band = [];
% test for psnr
i = 0;
for t = T
    i = i+1;
    QTt = QT(:,:,i);
    Thetat = Theta(:,:,i);
    % direct transform
    disp('Computing forward bandelet transform.');
    [MB,Rbandg] = perform_wavelet_bandelet_transform(M,Jmin,QTt,Thetat,1, options);
    % quantize
    MBt = perform_quantization(MB,t);
    MWt = perform_quantization(MW,t);
    
    % evaluate the nbr of bits via entropy
    Rwav  = evaluate_nbr_bits_wavelets(MW, Jmin, t);
    Rband = evaluate_nbr_bits_wavelets(MB, Jmin, t);

    disp('Computing inverse bandelet transform.');
    [Mb,Rbandg] = perform_wavelet_bandelet_transform(MBt,Jmin,QTt,Thetat,-1, options);
    
    disp('Computing inverse bandelet transform.');
    Mw = perform_wavelet_transform(MWt,Jmin,-1);

    psnr_wav = [psnr_wav, psnr( Mw,M )];
    psnr_band = [psnr_band, psnr( Mb,M )];
    bit_wav = [bit_wav, Rwav];
    bit_band = [bit_band, Rband+Rbandg];
end

bit_wav = bit_wav/n^2;
bit_band = bit_band/n^2;
plot(bit_wav, psnr_wav, bit_band, psnr_band);
axis tight;
legend('Wavelets', 'Bandelets');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m1 = min( min(bit_wav), min(bit_band) );
m2 = max( max(bit_wav), max(bit_band) );
bpp = linspace(m1,m2,200);

interp_psnr_wav = interp1(bit_wav, psnr_wav, bpp);
interp_psnr_band = interp1(bit_band, psnr_band, bpp);

interp_gain = interp_psnr_band - interp_psnr_wav;
I = find(~isnan(interp_gain));
interp_gain = interp_gain(I);
bpp = bpp(I);
plot(bpp, interp_gain);
title('PSNR Gain wav->band');
xlabel('bbp');
ylabel('+PSNR');
axis tight;

if save_images
    if wavelet_vm==4
        str_wav = '7-9';
    elseif wavelet_vm==0
        str_wav = 'haar';
    end
    if options.use_single_qt==1
        str_geom = 'singlegeom';
    else
        str_geom = 'triplegeom';
    end
    str = [name '_psnrgain_' str_wav '_' str_geom];
    saveas(gcf, [str '.png'], 'png');
end