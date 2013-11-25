function [quality,ImageMean,ImageMax] = SCIELAB_FR(origImg, distImg)

% origImg = imread('hats.bmp');
% distImg = imread('hats_JPEG.bmp');

%% 2.  Load the calibration information
% SCIELAB calculates the difference between two images, specified in lms or xyz units
% Here we use the spectral power distribution of the red, green and blue color primaries in a display (displaySPD)
% and the spectral sensitivities of the cone photodetectors (spCones) 
% to calculate a transform that maps linear display rgb values into lms units 

path(path, 'D:\tools\matlab\workspace\Image Quality Metrics\SCIELAB_FR');

% we do this over a range of wavelengths ranging between 400 and 700 nn,
% sampled in steps of 10 nm
wave = 400:10:700;

% addpath /Image Quality Metrics/SCIELAB
% read in the spectral sensitivities of the cone photodetectors
coneFile = fullfile('SCIELAB_FR','SmithPokornyCones');
spCones = vcReadSpectra(coneFile,wave);   %plot(wave,spCones)

% read in the spectral power of the red, green and blue color channels of an sRGB display
displayFile = fullfile('SCIELAB_FR','crtSPD');
displaySPD = vcReadSpectra(displayFile,wave);   %plot(wave,displaySPD)

% transform from linear display rgb to lms units
rgb2lms = spCones'* displaySPD;

% The display gamma will be used to convert display DAC rgb values into
% LINEAR display rgb values
displayGamma = load('displayGamma');

% SCIELAB requires that you define a white point 
% this is the lms or xyz values for what you define as white.
% here we assume that the maximum linear rgb values (1 1 1) define the
% white point
% note that display DAC rgb values range between 0 and 255 and 
% the LINEAR display rgb values range between 0 and 1
rgbWhite = [1 1 1];
whitePt = rgbWhite * rgb2lms';

%% 3.1  -- Convert the RGB data to LMS (or XYZ if you like).
% SCIELAB accepts either lms or xyz values .. here we use lms values

imgRGB  = dac2rgb(origImg,displayGamma.gamma);
img1LMS = changeColorSpace(imgRGB,rgb2lms);

imgRGB = dac2rgb(distImg,displayGamma.gamma);
img2LMS = changeColorSpace(imgRGB,rgb2lms);

imageformat = 'lms';

%% 4. --  Run the scielab function.

% This parameter determines how many image samples per degree
% This can depend on the viewing distance.  So we set the expected visual
% angle of the image, and we then compute how many spatial samples per degree.
horizontalAngle = 5;  % degrees  
sz = size(origImg); % number of pixels in 5 degrees
sampPerDeg = round(sz(2)/horizontalAngle); % samples per degree of visual angle
% Note: The horizontalAngle is the visual angle subtended by the image at a specified viewing distance 
%       Visual angle in units of radians is 2 arctan(S/2D) where S is the size of the image and D is the viewing distance
%       degrees = radians (180/3.1416)
%       so visual angle in units of degrees is 2 arctan(S/2D) * (180/3.1416)
%       (e.g. size of image on the display is 0.5 cm and viewing distance is 5.71 cm, then visual angle in degrees is ~5

% Run CIELAB 2000
params.deltaEversion = '2000';
params.sampPerDeg = sampPerDeg;
params.imageFormat = imageformat;
params.filterSize = sampPerDeg;
params.filters = [];
% params.filterversion = 'original';
params.filterversion = 'distribution';
errorImage = scielab(img1LMS, img2LMS, whitePt, params);
ImageMean = mean(mean(errorImage));
ImageMax = max(max(errorImage));

quality = mean(mean(errorImage));
% [m,n] = size(errorImage);
% D0=0.8;
% %%%%%
% k =sum(sum(abs(errorImage)))/m*n;
% % quality =sum(sum(abs(errorImage)))/m*n;
% % load k_contourlet;
% % k_contourlet = [k_contourlet k];
% % save k_contourlet k_contourlet;
% % k =sum(1-((NC2-NC1)./NC1).^2);
% % quality = D1/(D1+log2(1+k/D0));
% % 07/20
% quality=1/(1+log2(1+k'/D0));


%% 5. --  Examining and interpreting the results.
% figure(1); clf
% imagesc(errorImage);  colorbar('vert');
% title('S-CIELAB error map')
% 
% % We think this is 200 for certain angles (near 10) of the image
% figure(2);
% hist(errorImage(:),50);
% fprintf('Number delta (E > 10): %.0f\n',sum(errorImage(:) > 10));   
% 
% % Look at the spatial distribution of the errors.
% gLevels = 128;
% errorTruncated = min(gLevels*(errorImage/10),gLevels*ones(size(errorImage)));
% 
% figure(3);
% image(errorTruncated); axis image; colormap(gray(gLevels));
