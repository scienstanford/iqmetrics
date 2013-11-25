%% s_scielabHarmonicsUncorrelatedError.m
%
% The purpose of this experiment is to simulate the stimulus conditions
% used in Foley and Brown, 1980.  We will see if SCIELAB can predict the
% contrast threholds.
%
%
% (c) Imageval 2012

%% Background
% Foley and Legge, 1980 screen subtended 10 degrees horizontally and 6 degrees vertically
% 2 AFC method - which interval contained the target
% for 2 cpd the 3 observers' threshold contrasts varied between 0.28 and 0.51
% Masker frequencies 1.0, 1.4, 1.7, 2.0, 2.4, 2.8, and 4.0 cpd.
% For each masking frequency, the 11 masking contrasts (percent contrasts) were:
% 0.05, 0.10, 0.20, 0.40,0.80, 1.6, 3.2, 6.4, 12.8, 25.6, and 51.2%.

% Figure 1: Results for narrow field condition: masker frequency was 2.8 cpd and the gratings were 0.75°
% wide
% narrow field masking: stimuli were 0.75 deg, again centered about a fixation point
% plot threholds (0.1 to 10 % contrast) on y axis and
% log masking contrast (0.1 to 300 %contrast) on the x axis

% Figure 2: Results for wide field condition
% wide field masking: stimuli were 6 degrees wide symmetric about a fixation point
% signal thresholds as a function of the contrast of maskers.
% Plot threholds (0.1 to 10 % contrast) on y axis and
% log masking contrast (0.1 to 50 % contrast) on x axis

%% Parameters of the target and mask
% Figure 1 shows threshold contrast as a function of masking contrast
% Mask was 2.8 and target was 2.0 cpd
% 0.75 degree field
%
% when MaskFreq = TargetFreq, this is a discrimination task
threshold = 5.0;
TargetFreq = 2.0; TargetContrast = 10/100;
MaskFreq = 2.0;
MaskContrast = [10/100, 20/100, 30/100, 40/100, 50/100, 60/100, 70/100, 80/100, 90/100];
% units are percent ... what are the contrast units in sceneCreate
sParams = scParams;
sParams.sampPerDeg = 128; % assuming FOV = 1 degree (narrow field of view was 0.75 deg)
meanL = 200; %
Params = imageHparams; % create default harmonic parameters
Params.ph  = 0;
Params.ang = 0;
Params.freq = TargetFreq;

%% Target
Params.contrast =  TargetContrast;
% Params.ang = pi/6;
% Params.GaborFlag = 0.2;
Target = sceneCreate('harmonic',Params);
Target = sceneSet(Target,'name','Target');
Target  = sceneAdjustLuminance(Target,meanL);
vcAddAndSelectObject(Target);
sceneWindow;
%% MaskPlusTarget
Params.freq(2) = MaskFreq;

for ii = 1: length(MaskContrast)
    Params.contrast(2) = MaskContrast(ii);
    Params.ph(2)  = 0;
    Params.ang(2) = 0;
    MaskPlusTarget = sceneCreate('harmonic',Params);
    MaskPlusTarget = sceneSet(MaskPlusTarget,'fov',1);
    MaskPlusTarget = sceneSet(MaskPlusTarget,'name','MaskPlusTarget');
    MaskPlusTarget = sceneAdjustLuminance(MaskPlusTarget,meanL);
    vcAddAndSelectObject(MaskPlusTarget);
    sceneWindow
    
    %% Filter the Target and TargetPlusMask with chromatic-spatial filters for opponent channels
    %   First we have to convert images to XYZ values
    %   then we send in the XYZ images along with the XYZ for the white point
    % ****** What should the white point be? ****
    TargetXYZ  = sceneGet(Target,'xyz');
    max(TargetXYZ(2,:))
    MaskPlusTargetXYZ = sceneGet(MaskPlusTarget,'xyz');
    max(MaskPlusTargetXYZ(2,:))
    WhiteXYZ = sceneGet(Target,'illuminantxyz'); %
    
    %%
    scP = scParams;
    N = 50;                    % Try values from 20 to 200, for example
    scP.sampPerDeg    = N;
    scP.filterSize    = N;
    % Image size
    fprintf('Image size (deg): %.3f\n',size(TargetXYZ,1)/scP.sampPerDeg);
    
    % These are unblurred renderings of the opponent colors images
    imgOpp = imageLinearTransform(TargetXYZ, colorTransformMatrix('xyz2opp', 10));
    imgNOpp = imageLinearTransform(MaskPlusTargetXYZ, colorTransformMatrix('xyz2opp', 10));
    
    % Now we prepare the filters and blur them
    [scP.filters, scP.support] = scPrepareFilters(scP);
    [imgFilteredXYZ, imgFilteredOpp] = scOpponentFilter(TargetXYZ,scP);
    
    [imgFilteredNXYZ, imgFilteredNOpp] = scOpponentFilter(MaskPlusTargetXYZ,scP);
    
    %% Find the error xyz2 = xyz1 + E;
    [r,c,w] = size(imgFilteredNXYZ);
    totalE = imgFilteredNXYZ - imgFilteredXYZ;
    
    % Break the error into parts. This will become a function that takes XYZ or
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
    %
    
    % Find the part of the error that looks like the original.  One way is to
    % solve for the scale factor that matches all of XYZ terms
    %
    %   Solve for alpha that minimizes
    %     SSE(totalE - alpha*img)
    %
    method = 'ls';    % Least-squares
    alpha = metricsMaskedError(imgFilteredXYZ(:),totalE(:),method);
    maskedE = alpha*imgFilteredXYZ(:);
    unmaskedE = totalE(:) - alpha*imgFilteredXYZ(:);
    % unmaskedE = totalE(:) - alpha*imgFilteredXYZ(:);
    
    % Ooops.  Because this is a global orthogonalization, we still have little
    % local bits that look like the image.  We need to do this over small
    % blocks of the image?  Or gaussian regions?  And then use L4 norm for
    % summing across regions?
    
    % vcNewGraphWin; imagesc(reshape(totalE + 0.5,r,c,w));
    % vcNewGraphWin; imagesc(reshape(unmaskedE + 0.5,r,c,w));
    % vcNewGraphWin; imagesc(reshape(maskedE + 0.5,r,c,w));
    
    % Now create a new error image that differs from the original just by this
    % decorrelated (unmasked) error
    unmaskedXYZ2 = imgFilteredXYZ(:) + unmaskedE(:);
    maskedXYZ2 = imgFilteredXYZ(:) + maskedE(:);
    
    % Shape it back into the image format, and display it
    unmaskedXYZ2 = reshape(unmaskedXYZ2,r,c,w);
    maskedXYZ2 = reshape(maskedXYZ2,r,c,w);
    
    % f = vcNewGraphWin;
    % set(f,'name','Images after SCIELAB')
    % subplot(2,2,1); imagescRGB(xyz2srgb(imgFilteredXYZ)); title('Target')
    % subplot(2,2,2); imagescRGB(xyz2srgb(imgFilteredNXYZ));title('Target + Mask')
    % subplot(2,2,3); imagescRGB(xyz2srgb(unmaskedXYZ2));   title('Target + Uncorrelated Part of Mask')
    % subplot(2,2,4); imagescRGB(xyz2srgb(maskedXYZ2));   title('Target + Correlated Part of Mask')
    % % subplot(2,2,4); plot(totalE(:),unmaskedE(:),'.');
    % % subplot(2,2,4);
    % % xlabel('Total error'), ylabel('Unmasked error'); grid on
    % % We are interested in computing SCIELAB for xyz1 and unmaskedXYZ2.
    
    %% Calculate DE
    params.deltaEversion = '2000';   % Which CIELAB version
    deltaEImage1 = scComputeDifference(imgFilteredXYZ,imgFilteredNXYZ,WhiteXYZ,params.deltaEversion);
    deltaEImage2 = scComputeDifference(imgFilteredXYZ,unmaskedXYZ2,WhiteXYZ,params.deltaEversion);
    % vcNewGraphWin; imshow(deltaEImageScaled);
    if ( mean(deltaEImage1(:)) < threshold)
        TargetContrast1(ii) = TargetContrast + 0.1;
    end
    if ( mean(deltaEImage2(:)) < threshold)
        TargetContrast2(ii) = TargetContrast + 0.1;
    end
end

%
% f = vcNewGraphWin;
% set(f,'name','DeltaE')
% title('Delta E');
% subplot(2,2,1); imshow(iescale(deltaEImage1)); title('Target + Mask');
% subplot(2,2,2); hist(deltaEImage1(:));
% subplot(2,2,3); imshow(iescale(deltaEImage2)); title('Target + Uncorrelated Mask');
% subplot(2,2,4); hist(deltaEImage2(:));
%
% f = vcNewGraphWin;
% set(f,'name','UncorrelatedMaskOnly')
% title('UncorrelatedMaskOnly');
% unmaskedE = reshape(unmaskedE,r,c,w);
% imshow(iescale(unmaskedE));
% max(unmaskedE(:))

% end;


plot(MaskContrast,sCIELAB_DE,'r');
hold on;
plot(MaskContrast, sCIELAB_uncorrelatedDE,'g');
xlabel('DeltaE');
ylabel('Contrast');
