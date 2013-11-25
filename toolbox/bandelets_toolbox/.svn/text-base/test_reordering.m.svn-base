% test for pixel re-ordering and test of different 1D-wavelet scheme
%
%   Copyright (c) 2005 Gabriel Peyré

name = 'parabola';
n = 256;

% you can play on these two parameter to check 
% the consticency of the scheme with respect
% to curvature and size.
c = 0.2;    % curvature of the parabola
P = 16;     % size of the square

save_images = 0;

% repertory to ouput display
rep = 'images/';
if save_images
    warning off;
    mkdir(rep);
    warning on;
end
% svg string for image
base_str = [rep name '_P' num2str(P) '_c' num2str( round(10*c) ) ];

Jmin = 3;
Jmax = log2(n)-1;

if 1 || ~exist('MW')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % image loading
    if strcmp(name, 'parabola')
        % skip high details
        options.c = c;    % curvature of the parabola
        M_original = load_image( name, 1024, options );
        disp('Creating filetered image.');
        MW = perform_wavelet_transform(M_original,Jmax+1);
        M_original = MW(1:n,1:n);
    else
        M_original = load_image( name, n );
    end
    M_original = 256 * rescale(M_original);
    disp('Performing wavelet transform.');
    MW = perform_wavelet_transform(M_original,Jmin);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% region of interest
switch lower(name)
    case 'lena'
        roi_center = [0.8022, 0.6767];
        %     roi_center = [0.875,0.665];
    case 'peppers'
        roi_center = [0.821,0.306];     % line
        %    roi_center = [0.21,0.136];     % curve
    case 'barb'
        roi_center = [0.478,0.62];
        roi_center = [0.662,0.388];     % left pant
    otherwise
        roi_center = [0.5,0.5];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select a sub square
j = Jmax;   % scale
q = 1;      % orientation
[selx,sely] = compute_quadrant_selection(j,q);  
MWj = MW(selx,sely);
roi_centern = floor(roi_center*(2^j-1))+1;
M = MWj( roi_centern(1)-P/2+1:roi_centern(1)+P/2, roi_centern(2)-P/2+1:roi_centern(2)+P/2 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find best geometry
T = 10;
s = 12;
disp('Computing geometry.');
[MB,theta] = compute_best_direction(M,T,s);
disp('Performing displays.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the geometry
x = -1:2/(P-1):1;
a = cos(theta)*[1,-1]*1.3;
b = sin(theta)*[1,-1]*1.3;
clf;
hold on;
imagesc(x,x,M);
line(b,a);
hold off;
colormap gray(256);
axis image; axis off;
if save_images
    saveas(gcf, [base_str, '_geometry.png'], 'png');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% warp points, plot curves
[t,u,v,ts,us,vs] = perform_warping(M,theta,P);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot 1D functions
clf;
subplot( 2,1, 1 );
plot( ts, vs, '.-' );
title('Real positions plot.');
axis tight;
subplot( 2,1, 2 );
plot( 1:length(vs), vs, '.-' );
title('Integer-rounded positions plot.');
axis tight;
if save_images
    saveas(gcf, [base_str, '_1D_functions.png'], 'png');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% approximation using various schemes
w = [];
e = [];
lgd = {};

% HAAR
y = perform_haar_transform(v);
w = [w, y]; e = [e, l2error(y)];
lgd{end+1} = 'Haar';

% 7-9
y = perform_wavelet_transform(v);
w = [w, y]; e = [e, l2error(y)];
lgd{end+1} = 'Wavelet 7-9';


% alpert
alpert_vm = [ [0;0], [1;0], [0;1], [1;1], [2;0], [0;2] ];
pos = [t,u]';
y = perform_alpert_transform_2d(v,pos,alpert_vm, 1);
w = [w, y]; e = [e, l2error(y)];
lgd{end+1} = 'Alpert 1,X,Y';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot approximation rates
sel = 1:(P^2)/2;
clf;
plot( log2(sel), log2(e(sel,:)) );
ylabel('log2(|f-f_M|^2)');
xlabel('M');
legend(lgd);
axis([0, log2(max(sel)), -1, log2(max(e(:)))]);
if save_images
    saveas(gcf, [base_str, '_approximation.png'], 'png');
end

if exist('ARTICLE_save_square_size')
    ARTICLE_save_square_size;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPS exports for the article
v1 = vs( max(end/2-64+1,1):min(end/2+64,length(vs)) );
clf;
stem(v1, 'filled' );
axis tight;
axis off;
if save_images
    saveas(gcf, [base_str, '_stem1.eps'], 'eps');
end

y = perform_wavelet_transform(v1, 0);
clf;
stem(abs(y), 'filled');
axis tight;
axis off;
if save_images
    saveas(gcf, [base_str,  'stem2.eps'], 'eps');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% color display of the image (sexy display for cover ;-)
q = 16;
M = M(end/2-q/2+1:end/2+q/2, end/2-q/2+1:end/2+q/2);
eta = 1/5;

Mr = M*0+1;
Mv = M*0+1;
Mb = M*0+1;

I = find(M<0);
a = rescale( -M(I).*(M(I)<0) );
a = a.^eta;
Mv(I) = 1-a; Mb(I) = 1-a;

I = find(M>0);
a = rescale( M(I).*(M(I)>0) );
a = a.^eta;
Mr(I) = 1-a; Mv(I) = 1-a;

Mi = ones(q,q,3);
Mi(:,:,1) = Mr;  Mi(:,:,2) = Mv; Mi(:,:,3) = Mb;

% does this image remind you of something ?
imagesc(Mi);
if save_images
    imwrite(Mi, [base_str, '_color.png'], 'png');
end