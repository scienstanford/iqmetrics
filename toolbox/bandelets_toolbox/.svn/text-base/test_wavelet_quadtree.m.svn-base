% test the quadtree construction in wavelet domain
%
%   Copyright (c) 2005 Gabriel Peyr?

% set to 0 if you don't want to save the images
save_images = 0;

rep = 'images/';
warning off;
mkdir(rep);
warning on;

n = 512;

options.sigma = 3;
M = load_image('lena', n);

if save_images
    warning off;
    imwrite( M/256, [rep 'original.png'], 'png' )
    warning on;
end

Jmin = 3;
Jmax = log2(n)-1;
j_min = 2;
j_max = 4;
T = 100;
disp('Computing quadtree.');
[QT,Theta] = compute_wavelet_quadtree(M,Jmin,T,j_min,j_max);


MW = perform_wavelet_transform(M,Jmin, 1);

% save QT and Theta
if save_images
    warning off;
    imwrite( (QT-min(QT(:)))/(max(QT(:))-min(QT(:))), [rep sprintf('qt_bandelet_wavelet_%d_m%d.png', T,m)], 'png' );
    imwrite( Theta/pi, [rep sprintf('theta_bandelet_wavelet_%d_m%d.png', T,m)], 'png' );
    warning on;

    
    str = [rep sprintf('quadtree_wavelet_bandelet_%d', T)];
    clf;
    plot_wavelet_quadtree(QT,Theta,Jmax,1,MW);
    saveas(gcf, str, 'png');

    clf;
    plot_wavelet_quadtree(QT,Theta,Jmax,1,MW,2);
    saveas(gcf, [str, '_geom'], 'eps');
    
    clf;
    plot_wavelet_quadtree(QT,Theta,Jmax,1,MW,3);
    saveas(gcf, [str, '_struct'], 'eps');
end

% direct transform
disp('Computing forward bandelet transform.');
[MB,m_geom] = perform_wavelet_bandelet_transform(M,Jmin,QT,Theta,1);

disp('Computing inverse bandelet transform.');
[M1,m_geom] = perform_wavelet_bandelet_transform(MB,Jmin,QT,Theta,-1);

disp( ['|f-B^-1(B(f)|^2 (should be 0 for a correct reconstruction) --> ' sprintf( '%f', sum( (M(:)-M1(:)).^2 ) )] );

% performing truncated reconstruction
MBt = MB .* (abs(MB)>T);	% thresholded transform
[M1,m_geom] = perform_wavelet_bandelet_transform(MBt,Jmin,QT,Theta,-1);
clf;
imagesc(M1);
% number of coefficients used
m = round( sum( abs(MB(:))>T )+m_geom );


% save image
if save_images
    warning off;
    imwrite( M1/256, [rep sprintf('bandelet_wavelet_%d_m%d.png', T,m)], 'png' );
    disp( sprintf('T=%d, #coefficients=%d', T, m) );
    warning on;
end