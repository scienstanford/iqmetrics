%%Program Description
%This program is the main entry of the application.
%This program calculates the difference between Image and video Quality metrics
% 04.10.2011

%% Clear Memory & Command Window
clc;
clear all;
close all;

%% Read Original & Distorted Images
origImg = imread('hats.bmp');
distImg = imread('hats_JPEG.bmp');

%%
% origImg1 = origImg(:,:,1); %R
% origImg2 = origImg(:,:,2); %G
% origImg3 = origImg(:,:,3); %B

%%  Initialization: If the input image is rgb, convert it to gray image
noOfDim = ndims(origImg);
if(noOfDim == 3)
    origImg = rgb2gray(origImg);
end

noOfDim = ndims(distImg);
if(noOfDim == 3)
    distImg = rgb2gray(distImg);
end


%% Size Validation
origSiz = size(origImg);
distSiz = size(distImg);
sizErr = isequal(origSiz, distSiz);
if(sizErr == 0)
    disp('Error: Original Image & Distorted Image should be of same dimensions');
    return;
end

% 
% Full-reference Image Quality Metrics

disp('Full-reference image quality Metrics');
%% [1] Mean Square Error(MSE)
MSE = MeanSquareError(origImg, distImg);
disp(' [1] Mean Square Error = ');
disp(MSE);
% 
%% [2]Peak Signal to Noise Ratio (PSNR)
PSNR = PeakSignaltoNoiseRatio(origImg, distImg);
disp('[2] Peak Signal to Noise Ratio = ');
disp(PSNR);
% % 
%% [3] Normalized Cross-Correlation (NCC)
NK = NormalizedCrossCorrelation(origImg, distImg);
disp('[3] MNormalized Cross-Correlation  = ');
disp(NK);
% % 
%% [4] Average Difference (AD)
AD = AverageDifference(origImg, distImg);
disp('[4] Average Difference  = ');
disp(AD);
% % 
%% [5] Structural Content (SC)
SC = StructuralContent(origImg, distImg);
disp('[5] Structural Content  = ');
disp(SC);
% % 
%% [6] Maximum Difference (MD)
MD = MaximumDifference(origImg, distImg);
disp('[6] Maximum Difference = ');
disp(MD);
% 
%% [7] Normalized Absolute Error (NAE)
NAE = NormalizedAbsoluteError(origImg, distImg);
disp('[7] Normalized Absolute Error = ');
disp(NAE);
% 
%% [8] Universal Image Quality Index (UQI)
[quality, quality_map] = UQI_FR(origImg, distImg);
disp('[8] Universal Image Quality Index = ');
disp(quality);
figure, imshow(quality_map),title('Universal Image Quality Index Quality Map');
% 
%% [9] Structral SIMilarity (SSIM)
[mssim, ssim_map] = SSIM_FR(origImg, distImg);
disp('[9] Structral SIMilarity = ');
disp(mssim);
figure, imshow(ssim_map),title('Structral SIMilarity Quality Map');
% 
%% [10] Multi-scale Structural Similarity Index (MS-SSIM)
mssim = MS_SSIM_FR(origImg,distImg);
disp('[10] Multi-scale Structural Similarity Index (MS-SSIM) = ');
disp(mssim);
%
% [11] SVD-based Image Quality Measure
[scaMeasure, graMeasure] = SVD_FR(origImg, distImg);
disp('[11] SVD-based Image Quality Measure = ');
disp(scaMeasure);
figure, imshow(graMeasure),title('SVD-based Image Quality Map');
%
% [12] Visual signal-to-noise ratio for digital images 
[res] = VSNR_FR(origImg, distImg);
disp('[12] Visual signal-to-noise ratio = ');
disp(res);
%
%% [13] Information Fidelity Criterion (IFC)
ifc = ifc_FR(origImg,distImg);
disp('[13] Information Fidelity Criterion  = ');
disp(ifc);
% % 
%% [14] Visual Information Fidelity (VIF) measure 
vif= VIF_FR(origImg,distImg);
disp('[14] Visual Information Fidelity (VIF) measure  = ');
disp(vif);

% [15] VIF based on pixel-wise
vifp =vifp_mscale(origImg,distImg);
disp('[15] pixel domain of Visual Information Fidelity (VIF) measure  = ');
disp(vifp);

% [16] DCTune_FR
Q = DCTune_FR(origImg,distImg);
disp(' [16] DCTune = ');
disp(Q);
% %
% [17] pixel-based JND (Just-Noticeable Difference) model
JND=JND_pixel();
disp(' [17] pixel-based JND (Just-Noticeable Difference) model = ');
disp(JND);
figure, imshow(JND),title('JND Quality Map') 
% %  
%% [18] DCT-based JND (Just-Noticeable Difference) model
JND=JND_dct(distImg);
disp(' [18] DCT-based JND (Just-Noticeable Difference) model = ');
disp(JND);
figure, imshow(JND),title('JND Quality Map') 
% %   
%% [19] Nonlinear Weighted Signal to noise ratio for additive noise. 
NQM = NQM_FR(origImg,distImg);
disp(' [19] Nonlinear Weighted Signal to noise ratio = ');
disp(NQM);
% %  
%% [20] content based metric
FI = CBM_FR(origImg,distImg);
disp(' [20] content based metric = ');
disp(FI);
%
%% [21] Information content weighted structural similarity measure (IW-SSIM)
%% [22] Information content weighted Mean Square Error (IW-MSE)
%% [23] Information content weighted Peak Signal to Noise Ratio (IW-PSNR)
   [iwssim iwmse iwpsnr]= IW_SSIM_MSE_PSNR_FR(origImg, distImg)
    disp(' [21]iwssim = ');
    disp(iwssim);
    disp(' [22]iwmse = ');
    disp(iwmse);
    disp(' [23]iwpsnr = ');
    disp(iwpsnr);
% %%
% [24] FeatureSIM
[FSIM, FSIMc] = FSIM_FR(origImg,distImg)
disp(' [24] FSIM_FR = ');
disp(FSIM);

% %% Reduced Reference Image Quality Metric
% disp('reduced reference image quality metric');
% 
% % [1] Reduced reference image quality assessment using a wavelet-domain natural image statistic model
%      distortion = WNISM_RR(origImg,distImg);
%      disp(' [1]WNISM = ');
%      disp(distortion);
%      
% %% [2] Image quality assessment based on wavelet domain
%     quality = Wavelet_RR(origImg,distImg);
%     disp(' [2] Wavelet_RR = ');
%     disp(quality);
% %     
% % % 
%% [3] Image quality assessment based on contourlet domain
% quality = Contourlet_RR(origImg,distImg);
% disp(' [3] Contourlet_RR = ');
% disp(quality);
%     
% % % [4] Image quality assessment based on WBCT domain
% quality = WBCT_RR(origImg,distImg);
% disp(' [4] WBCT_RR = ');
% disp(quality);
% % 

% 
% %% No reference image quality metric
% disp('No reference image quality metric');
% 
% %% [1] No-Reference Perceptual Quality Assessment of JPEG Compressed Images
% score = JPEG_S_NR(distImg);
% disp(' [1] content based metric = ');
% disp(score);
% 
% %% [2] No-Reference Image Quality Assessment based on Discrete Cosine Transformation
% proportion = DCT_NR(distImg);
% disp(' [2] DCT_NR = ');
% disp(proportion);
% % 
% 
% %% [4] A Modular Framework for Constructing Blind Universal Quality Indices
% % [quality probs] = BIQI_NR(distImg);
% % disp(' [4] BIQI_NR = ');
% % disp(quality);














