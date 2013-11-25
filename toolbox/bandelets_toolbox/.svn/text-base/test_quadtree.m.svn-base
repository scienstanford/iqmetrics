% test the quadtree construction
%
%   Copyright (c) 2005 Gabriel Peyré


% set to 0 if you don't want to save the images
save_images = 0;

rep = 'images/';
warning off;
mkdir(rep);
warning on;

% size of the image
n = 128;
% number of VM of the transforms
global wavelet_vm;
wavelet_vm = 0;

name = '3contours';
M = load_image(name, n);

j_min = 2;
j_max = 4;
T = [1,10,100];
disp('Computing quadtree.');
[QT,Theta] = compute_quadtree(M,T,j_min,j_max);

% save QT and Theta
if save_images
    warning off;
    imwrite( (QT-min(QT(:)))/(max(QT(:))-min(QT(:))), [rep sprintf('qt_bandelet_space_%d_m%d.png', T)], 'png' );
    imwrite( Theta/pi, [rep sprintf('theta_bandelet_space_%d%d.png', T)], 'png' );
    warning on;

    str = [rep sprintf('quadtree_bandelet_space_%d', T)];
    clf;
    plot_quadtree(QT,Theta, M);
    saveas(gcf, str, 'png');

    clf;
    plot_quadtree(QT,Theta, M, 2);
    saveas(gcf, [str, '_geom'], 'eps');
    
    clf;
    plot_quadtree(QT,Theta, M, 3);
    saveas(gcf, [str, '_struct'], 'eps');
end

% choose one of the thresholds
k = 3; T = T(k);

% direct transform
disp('Computing forward bandelet transform.');
[MB,m_geom] = perform_bandelet_transform(M,QT(:,:,k),Theta(:,:,k),1);
disp( ['|f|^2 - |B(f)|^2 (should be 0 for an orthogonal transform) --> ' sprintf( '%f', sum(M(:).^2 - MB(:).^2) )] );

disp('Computing inverse bandelet transform.');
[M1,m_geom] = perform_bandelet_transform(MB,QT(:,:,k),Theta(:,:,k),-1);

disp( ['|f-B^-1(B(f)|^2 (should be 0 for a correct reconstruction) --> ' sprintf( '%f', sum( (M(:)-M1(:)).^2 ) )] );

% performing truncated reconstruction
MBt = MB .* (abs(MB)>T);	% thresholded transform
[M1,m_geom] = perform_bandelet_transform(MBt,QT(:,:,k),Theta(:,:,k),-1);
clf;
imagesc(M1);
% number of coefficients used
m = round( sum( abs(MB(:))>T )+m_geom );


% save image
if save_images
    warning off;
    imwrite( M1/256, [rep sprintf('bandelet_space_%d_m%d.png', T,m)], 'png' );
    disp( sprintf('T=%d, #coefficients=%d', T, m) );
    warning on;
end