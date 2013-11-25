% s_imageCorrelatedError
%
% We are dividing image error into two parts - one that is a scaled version
% of the image and a second that is orthogonal to the image.
% 
% We think that each error term should be considered as having different
% weights - the scaled version of the original is masked and the orthogonal
% error is handled by S-CIELAB.
%
% In this script we are doing some experiments to test the idea and
% illustrate it.
%
% Question: Because this is a global orthogonalization, we still have little
% local bits that look like the image.  Should we calculate correlations for small local regions
% in the image?  

% (c) Imageval, 2011

%% Read an image, and distort it, separate out the two error parts
%
% Show some images of this and consider how the correlated image error
% might be reduced because of adaptation, while the uncorrelated image has
% a larger perceptual error per root mean squared errror value because,
% well, it is not masked.

%% Clear variables
s_initISET


%% Start with simple RGB values, nothing calibrated, as an illustration
fname    = 'eagle.jpg';
% fname    = 'hImage.jpg';
% fname    = 'mcc3.jpg';
% fname    = 'hats.jpg';
fullName = fullfile(isetRootPath,'data','images','rgb',fname);
img      = imread(fullName);
% imshow(img)

%% Create a noisy version of the image


% dMethod = 'Scale Contrast';
% dMethod = 'Gaussian Noise';
dMethod = 'JPEG Compress';
switch(dMethod)
    case 'Scale Contrast'
        img  = double(img);
        imgN = img;
        imgN = imageDistort(imgN,dMethod,0.1);
        mx1 = max(img(:));
        mx2 = max(imgN(:));
        mx  = max(mx1,mx2);
    case 'Gaussian Noise'
        img  = double(img);
        imgN = img;
        imgN = imageDistort(imgN,dMethod,50);
        mx1 = max(img(:));
        mx2 = max(imgN(:));
        mx  = max(mx1,mx2);
    case 'JPEG Compress'
        imgN = imageDistort(img,dMethod,30);
        img = double(img);
        imgN = double(imgN);
        mx  = 255;
end

% vcNewGraphWin; imshow(img)
% vcNewGraphWin; imshow(imgN)

%% Calculate XYZ values

% Suppose the original image is an sRGB image (rgb values designed to be
% displayed on a linear(?) sRGB display)
imgXYZ  = srgb2xyz(img/mx);  
% vcNewGraphWin; imagescRGB(double(img));
% vcNewGraphWin; imagescRGB(xyz2srgb(imgXYZ));

imgNXYZ = srgb2xyz(imgN/mx);

whiteXYZ = srgb2xyz(ones(1,1,3)); % max(Y(:)); whiteXYZ(2)
% max(imgXYZ(:))
% max(imgNXYZ(:))
% whiteXYZ


%% Set up parameters to spatially filter using SCIELAB opponent channels

% We will experiment with different samples per degree
% Changing the samplesPerDeg effectively changes the FOV
% High quality might be 250, low quality 80

% When N = 40, the image spans 6.4 degrees
% It has no real high frequency  content (beyond 20 cpd)
% So, there isn't much to blur
% When N = 250, we are simulating a 1 deg image. In that case, the image
% content has spatial frequencies up to 125 c/deg and thus the blurred
% version is very different from the original.
scP = scParams;
N = 50;                    % Try values from 20 to 200, for example
scP.sampPerDeg    = N;   
scP.filterSize    = N;
% Image size
fprintf('Image size (deg): %.3f\n',size(imgXYZ,1)/scP.sampPerDeg)

%% Filter using SCIELAB opponent channels
% We first transform XYZ images into opponent channel images
% and then we blur - each channel has a different spatial filter
% and then transform from opponent space back to XYZ 

% These are unblurred renderings of the opponent colors images
gam = 0.4;
imgOpp = imageLinearTransform(imgXYZ, colorTransformMatrix('xyz2opp', 10));
% imagescOPP(imgOpp,gam);
imgNOpp = imageLinearTransform(imgNXYZ, colorTransformMatrix('xyz2opp', 10));

% Now we prepare the filters and blur them
[scP.filters, scP.support]       = scPrepareFilters(scP);
[imgFilteredXYZ, imgFilteredOpp] = scOpponentFilter(imgXYZ,scP);
imgFilteredXYZ = ClipXYZImage(imgFilteredXYZ,whiteXYZ);
% imagescOPP(imgFilteredOpp,gam);

[imgFilteredNXYZ, imgFilteredNOpp] = scOpponentFilter(imgNXYZ,scP);
imgFilteredNXYZ = ClipXYZImage(imgFilteredNXYZ,whiteXYZ);
% imagescOPP(imgFilteredNOpp,gam);

%% Find the total error in XYZ coordinates 

[r,c,w] = size(imgFilteredNXYZ);
totalE = imgFilteredNXYZ - imgFilteredXYZ;

%% Break the error into parts. 

% This will become a function that takes XYZ or
% OPP and the error and returns a scale factor for the correlated error.
% The scale factor could be computed various different ways.  We then
% calculate the uncorrelated (or unmasked) error as per below.
%
% Note that the error is computed AFTER the blurring from the spatial part
% of the calculation.  That is probably important, and we should check the
% significance of this with images.
%
% Create a new error term that removes the masked component from the total
% error.  We express
%  totalE = unmaskedE + alpha*img
%  unmaskedE = totalE - alpha*img

% Find the part of the error that looks like the original (masked error)  

% One way is to
% solve for the scale factor that matches all of XYZ terms
%
%   Solve for alpha that minimizes
%     SSE(totalE - alpha*img)
%
method = 'ls';    % Least-squares
Y = imgFilteredXYZ(:,:,2);
YE = totalE(:,:,2);
% This is the total corelation estimate
alpha = metricsMaskedError(imgFilteredXYZ(:),totalE(:),method); 

% We can calculate the correlation for luminance only
% Joyce thinks this is a bad idea, but we try it anyway.
% alpha = metricsMaskedError(Y(:),YE(:),method); 
% sprintf('Correlated error %f\n',alpha)

% The masked error is like the image
maskedE = alpha*imgFilteredXYZ(:);

% The unmasked error is the part of the error that is unlike the original
unmaskedE = totalE(:) - alpha*imgFilteredXYZ(:);
% unmaskedE = totalE(:) - alpha*imgFilteredXYZ(:);

% vcNewGraphWin; imagesc(reshape(totalE + 0.5,r,c,w));
% vcNewGraphWin; imagesc(reshape(unmaskedE + 0.5,r,c,w));
% vcNewGraphWin; imagesc(reshape(maskedE + 0.5,r,c,w));

%% Create  error images that differ  by only unmasked or masked error

% Unmasked only
unmaskedXYZ2 = imgFilteredXYZ(:) + unmaskedE(:);

% Masked only
maskedXYZ2 = imgFilteredXYZ(:) + maskedE(:);

% Shape it back into the image format, and display it
unmaskedXYZ2 = reshape(unmaskedXYZ2,r,c,w);
maskedXYZ2   = reshape(maskedXYZ2,r,c,w);

f = vcNewGraphWin; 
set(f,'name','Images after SCIELAB')
subplot(2,2,1); imagescRGB(xyz2srgb(imgFilteredXYZ)); title('Image')
subplot(2,2,2); imagescRGB(xyz2srgb(imgFilteredNXYZ));title('Image + all error ')
subplot(2,2,3); imagescRGB(xyz2srgb(unmaskedXYZ2));   title('Image + only uncorrelated error ')
subplot(2,2,4); imagescRGB(xyz2srgb(maskedXYZ2));     title('Image + only correlated error ')
% subplot(2,2,4); plot(totalE(:),unmaskedE(:),'.'); 
% subplot(2,2,4); 
% xlabel('Total error'), ylabel('Unmasked error'); grid on
% We are interested in computing SCIELAB for xyz1 and unmaskedXYZ2.

%% Calculate S-CIELAB DE

params.deltaEversion = '2000';   % Which CIELAB version
deltaEImage1 = scComputeDifference(imgFilteredXYZ,imgFilteredNXYZ,whiteXYZ,params.deltaEversion);
deltaEImage2 = scComputeDifference(imgFilteredXYZ,unmaskedXYZ2,whiteXYZ,params.deltaEversion);
deltaEImage3 = scComputeDifference(imgFilteredXYZ,maskedXYZ2,whiteXYZ,params.deltaEversion);

%% Compare the dE values

f = vcNewGraphWin; 
set(f,'name','DeltaE')
title('Delta E');
rr = 3;
cc = 2;
subplot(rr,cc,1); imshow(ieScale(deltaEImage1)); title('DE for Image + Total error');
subplot(rr,cc,2); hist(deltaEImage1(:));

subplot(rr,cc,3); imshow(ieScale(deltaEImage2)); title('DE for Image + Uncorrelated error');
subplot(rr,cc,4); hist(deltaEImage2(:));

subplot(rr,cc,5); imshow(ieScale(deltaEImage3)); title('DE for Image + Correlated error');
subplot(rr,cc,6); hist(deltaEImage3(:));

%% Show the error images 

f = vcNewGraphWin; 
set(f,'name',' Error Comparison')
subplot(1,3,1)
totalE = reshape(maskedE + unmaskedE,r,c,w);
imshow(ieScale(totalE));
title('Total  ');

subplot(1,3,2)
unmaskedE = reshape(unmaskedE,r,c,w);
imshow(ieScale(unmaskedE));
title('Uncorrelated  ');
% max(unmaskedE(:))

subplot(1,3,3)
maskedE = reshape(maskedE,r,c,w);
imshow(ieScale(maskedE));
title(sprintf('Correlated  (alpha = %f)',alpha));

%%  Which error predicts the total error?  Masked or unmasked?

vcNewGraphWin;
subplot(1,2,1)
plot(deltaEImage1(:),deltaEImage2(:),'.')
xlabel('Total dE');
ylabel('Unmasked dE')

subplot(1,2,2)
plot(deltaEImage1(:),deltaEImage3(:),'.')
xlabel('Total dE');
ylabel('Masked dE')

