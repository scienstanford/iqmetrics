%%  s_chrominanceChangesOnly
%
%   Purpose: Distortions to images that do not affect luminance will appear
%   to be unchanged to metrics that only rely on changes in luminance (e.g.
%   SSIM).  This script illustrates this by keeping luminance constant but
%   changing other components.
%   Normally, changes in chrominance will not affect image quality so this
%   is not a problem. For example, we can blur the chrominance channels
%   quite a lot without noticing the differences. Hence, a metric that does
%   not consider distortions in the chrominance channels may still be
%   useful.
%   Although we may be able to demonstrate that distortions in chrominance
%   can affect image quality, the conditions in which we create these
%   distortions may not occur in normal image processing pipelines.  In
%   other words, it may be difficult to find cases where changes in an
%   image processing pipeline affect chrominance but not luminance.
%   Nonetheless, this script illustrates that changes in chrominance,
%   however unlikely, can affect image quality and metrics that rely on
%   luminance changes only will not be sensitive to these differences.

%       

%%  RGB2XYZ
RGB = imread('hats.jpg');
RGB1 = double(imread('hats.jpg'));

XYZ1 = rgb2xyz(RGB1); % we should replace this with an ISET function

%% Change X channel
X1 = XYZ1(:,:,1); 
% The Arnold function scrambles an image, we use this to scramble the
% pixels in the X channel
X3= arnold(X1,1);
[m,n] = size(X1);
% x1 = X1-20;  %  0.1~1
%X3 = X1+randn([m,n])*5; % the different x can influence the result of SSIM
Y1 = XYZ1(:,:,2); 
% Y3 = Y1/100;
%y1 = Y1+randn([m,n])*200;
%y1 = Y1-60;
Z1 = XYZ1(:,:,3);
% z1 = Z1+randn([m,n])*100; % <120
%  z1 = Z1+10;
XYZ2(:,:,1)=X3;
XYZ2(:,:,2)=Y1; 
XYZ2(:,:,3)=Z1;

plot(X1,X3,'k-');

RGB2 =xyz2rgb(XYZ2);
% note that some values are negative or > 255 indicating that the values are out of
% gamut  


%% Now convert to XYZ again to see that Y is unchanged
XYZ3 = rgb2xyz(RGB2);  % THIS IS NOT A CORRECT TRANSFORM

X2 = XYZ3(:,:,1); 
Y2 = XYZ3(:,:,2);
Z2 = XYZ3(:,:,3); 

figure(1);
plot(Y1,Y2,'r-');
xlabel('Y1');
ylabel('Y2');
figure(2);
plot(X1,X2,'c-');
xlabel('X1');
ylabel('X2');
figure(3);
plot(Z1,Z2,'k-');
xlabel('Z1');
ylabel('Z2');
hold on;

%% problem is that some values are now out of gamut
% How do I handle this ... it is a real problem for SCIELAB

RGB3 = ieClip(RGB2,0,1) * 255; 


%% Compare RGB1 and RGB2 using SSIM metric

% Initialization: If the input image is rgb, convert it to gray image
origImg = RGB1;
distImg = RGB2;
noOfDim = ndims(origImg);
if(noOfDim == 3)
    origImg = rgb2gray(origImg);
end

noOfDim = ndims(distImg);
if(noOfDim == 3)
    distImg = rgb2gray(distImg);
end



[mssim, ssim_map] = SSIM_FR(origImg, distImg);
disp('[9] Structral SIMilarity = ');
disp(mssim);
figure, imshow(ssim_map),title('Structral SIMilarity Quality Map');

% Note that the conversion of rgb2gray does introduce a little luminance
% differences

%% SCIELAB
% some important viewing conditions
vDist = 0.3;          % 12 inch viewing distance
dispCal = 'crt.mat';   % Calibrated display

% the scielab calculation (again, see s_scielabExample.m for a more detailed calculation
[eImage,s1,s2] = scielabRGB(RGB1, RGB2, dispCal, vDist);

% This is the mean delta E
mean(eImage(:))

% Show the RGB images as scenes, illustrating how the RGB data were
% converted to SPDs using the calibrated display
vcAddAndSelectObject(s1); vcAddAndSelectObject(s2);sceneWindow;

% Examining and interpreting the results.
vcNewGraphWin;
imagesc(eImage);
colorbar('vert');
title('S-CIELAB error map')

vcNewGraphWin;
hist(eImage(:),100)
title('S-CIELAB delta E histogram')

% calculate the mean DE for values greater than 2
count = 0;
DEdifs = 0;
rows = max(size(eImage(:,1)));
cols = max(size(eImage(1,:)));
for ii = 1:max(rows)
    for jj = 1: max(cols)
        if eImage(ii,jj) > 2.0
            DEdifs = eImage(ii,jj) + DEdifs;
             count = count +1;
    end
    end
end
MeanAbove2 = mean(DEdifs/count);
percent = (count/ (rows * cols ))* 100;
fprintf('Mean of Delta E with values greater than 2: %f\n',MeanAbove2);
fprintf('Percent of Delta E with values greater than 2: %f\n',percent);