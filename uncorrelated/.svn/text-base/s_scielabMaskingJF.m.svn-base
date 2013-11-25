%% s_scielabFreqMasking.m
%
% The purpose of this experiment is to characterize how S-CIELAB predicts
% the probability that people can detect a target superimposed on a
% background mask 
%
% The more a target looks like a mask, the more difficult it is to detect
% the target When the target has the same spatial frequency as the mask,
% the task is also called "contrast discrimination"
%
%  Set a target to be a particular frequency
%  Set the mask to be a particular frequency
%  Calculate DE for (target plus mask) vs. (mask alone)
%
%
% (c) Imageval 2012

%% Notes

% we are treating the stimuli as scenes ... should we treat them as displayed images?
% do we want to have a different structure in ISET that is the displayed image
% (i.e. display pixel structure, spd, gamma)

%% Clear variables
s_initISET

%% Set variables

% sCIELAB params:
scP = scParams;
scP.sampPerDeg = 128; % assuming FOV = 1 degree (narrow field of view was 0.75 deg)
meanL = 200;              % Mean luminance for a bright display

% Set harmonic parameters for two harmonics.
% We set the target contrast to zero for the mask alone case.
% We vary the target contrast for mask plus target to measure threshold.
% We use the first term as the mask and the 2nd as the target.
Params = imageHparams;
Params.ph(1)  = 0;
Params.ang(1) = 0;
Params.ph(2)  = 0;
Params.ang(2) = 0;
% Params.ang = pi/6;
% Params.GaborFlag = 0.2;

%% Initialize harmonic parameters
TargetFreq = 2.0;
MaskFreq = 2.0;
Params.freq(1) = MaskFreq;
Params.freq(2) = TargetFreq;

Params.contrast(1) = 0.0;   % Contrast units are [0,1]

%%
MaskContrast = (0:0.3: 0.9);
TargetContrast = (0:.3:0.9);
meanDE = zeros([length(TargetContrast), length(MaskContrast)]);

%% Create the mask alone stimulus 
for jj = 1:length(MaskContrast)
    Params.contrast(1) = MaskContrast(jj);
    Params.contrast(2) = 0;
    Mask  = sceneCreate('harmonic',Params);
    Mask  = sceneSet(Mask,'name','Mask');
    Mask  = sceneAdjustLuminance(Mask,meanL);
    Mask  = sceneSet(Mask,'fov',1);
    vcAddAndSelectObject(Mask);
    sceneWindow;
    
    %% Filter in color opponent space
    
    scP = scParams;
    scP.sampPerDeg    = sceneGet(Mask,'row');
    scP.filterSize    = sceneGet(Mask,'row');
    
    MaskXYZ = sceneGet(Mask,'xyz');
    % Go into color opponent space
    % imgOpp = imageLinearTransform(MaskXYZ, colorTransformMatrix('xyz2opp', 10));
    % And blur
    [scP.filters, scP.support] = scPrepareFilters(scP);
    MaskFilteredXYZ = scOpponentFilter(MaskXYZ,scP);
    
    % vcNewGraphWin; imagescRGB(xyz2srgb(MaskFilteredXYZ));
    
    %% Create the second stimulus and filter in color opponent space
    % Scene two is the target plus the mask, with higher target mask
    
    for ii = 1:length(TargetContrast)
        Params.contrast(2) = TargetContrast(ii);
        MaskTarget = sceneCreate('harmonic',Params);
        MaskTarget = sceneSet(MaskTarget,'name','MaskTarget');
        MaskTarget = sceneAdjustLuminance(MaskTarget,meanL);
        MaskTarget = sceneSet(MaskTarget,'fov',1);
        vcAddAndSelectObject(MaskTarget);
        sceneWindow;
        
        MaskTargetXYZ = sceneGet(MaskTarget,'xyz');
        MaskTargetFilteredXYZ = scOpponentFilter(MaskTargetXYZ,scP);
        
        % Calculate DE
        % What should WhiteXYZ be set to?
        %     max(MaskXYZ(2,:))
        %     max(MaskTargetXYZ(2,:))
        
        % We set the harmonic as if the mean is a gray scale (20%) card.
        WhiteXYZ = sceneGet(MaskTarget,'illuminantxyz'); %
        
        deltaEImage1 = scComputeDifference(MaskFilteredXYZ,MaskTargetFilteredXYZ,WhiteXYZ,scP.deltaEversion);
        % deltaEImage2 = scComputeDifference(MaskFilteredXYZ,unmaskedXYZ2,WhiteXYZ,scP);
        % vcNewGraphWin; imshow(deltaEImageScaled);
        
        meanDE(ii,jj) = mean(deltaEImage1(:));
        % vcNewGraphWin; imagesc(deltaEImage1); colormap(gray)
        % vcNewGraphWin; plot(deltaEImage1(64,:),'k-'); grid on
        %     f = plotScene(Mask,'hline luminance',[1 64]); uDataM = get(f,'userdata');
        %     g = plotScene(MaskTarget,'hline luminance',[1 64]); uDataMT = get(g,'userdata');
        %     vcNewGraphWin; plot(uDataM.pos,uDataM.data - uDataMT.data)
        
        % scielab2(ii,jj) = mean(deltaEImage2(:));
    end
end

%%

vcNewGraphWin
plot(TargetContrast,meanDE(:,1),'r*');
pause
hold on
plot(TargetContrast,meanDE(:,2),'b*');
pause
plot(TargetContrast,meanDE(:,3),'g*');
pause
plot(TargetContrast,meanDE(:,4),'b*');
xlabel('Target Contrast');
ylabel('DeltaE');



% plot(TargetContrast, scielab2,'g');

%%
% This should be a function that accepts Image1 and Image2 and returns two DE values,
% one based on basic SCIELAB and one based on UC_SCIELAB (uncorrelated)
% a switch to display
% Image1
% Image2
% Image1 plus the part of Image2 that is correlated with Image1
% Image1 plus the part of Image2 that is uncorrelated with Image1
% deltaE images:
%   'DE for Image + Total error'
%   'DE for Image + Uncorrelated error'
%   'DE for Image + Correlated error'
% There is the issue of what WhiteXYZ should be

% Find the error xyz2 = xyz1 + E;
%     [r,c,w] = size(MaskFilteredXYZ);
%     totalE = MaskTargetFilteredXYZ - MaskFilteredXYZ;
%     method = 'ls';    % Least-squares
%     alpha = metricsMaskedError(MaskFilteredXYZ(:),totalE(:),method);
%     maskedE = alpha*MaskFilteredXYZ(:);
%     unmaskedE = totalE(:) - alpha*MaskFilteredXYZ(:);

% Now create a new error image that differs from the original just by this
% decorrelated (unmasked) error
%     unmaskedXYZ2 = MaskFilteredXYZ(:) + unmaskedE(:);
%     maskedXYZ2 = MaskFilteredXYZ(:) + maskedE(:);
%
%     % Shape it back into the image format, and display it
%     unmaskedXYZ2 = reshape(unmaskedXYZ2,r,c,w);
%     maskedXYZ2 = reshape(maskedXYZ2,r,c,w);

%% Keep target frequency, mask frequency and target contrast constant - vary mask contrast


%%                                                              %%%%%%
%% Keep target frequency, target contrast, mask contrast constant, vary mask frequency


%% Clear variables
% s_initISET
% % sCIELAB params:
% scP = scParams;
% scP.sampPerDeg = 128; % assuming FOV = 1 degree (narrow field of view was 0.75 deg)
% meanL = 200; %
% % Set frequency parameters
% Params = imageHparams; % create default harmonic parameters
% Params.ph(1)  = 0;
% Params.ang(1) = 0;
% Params.ph(2)  = 0;
% Params.ang(2) = 0;
% % Params.ang = pi/6;
% % Params.GaborFlag = 0.2;
% 
% % when MaskFreq = TargetFreq, this is a discrimination task
% % units are percent ... what are the contrast units in sceneCreate
% %% Keep target frequency, mask frequency and mask contrast constant - vary target contrast
% % Scene one is the target plus mask at one contrast
% TargetFreq = 2.0; 
% MaskFreq = 2.0;  % when MaskFreq = TargetFreq, we should observer Weber's law
% Params.freq(1) = TargetFreq;
% Params.freq(2) = MaskFreq;
% Params.contrast(1) = 0.4;
% Params.contrast(2) = 0.4;
% Mask = sceneCreate('harmonic',Params);
% Mask = sceneSet(Mask,'name','Mask');
% Mask  = sceneAdjustLuminance(Mask,meanL);
% vcAddAndSelectObject(Mask);
% sceneWindow;
% scP = scParams;
% N = 50;                    % Try values from 20 to 200, for example
% scP.sampPerDeg    = N;
% scP.filterSize    = N;
% % Image size
% fprintf('Image size (deg): %.3f\n',size(Mask,1)/scP.sampPerDeg);
% 
% MaskXYZ = sceneGet(Mask,'xyz');
% % Go into color opponent space
%  imgOpp = imageLinearTransform(MaskXYZ, colorTransformMatrix('xyz2opp', 10));
%  % And blur
%  [scP.filters, scP.support] = scPrepareFilters(scP);
%      [MaskFilteredXYZ, MaskFilteredOpp] = scOpponentFilter(MaskXYZ,scP);
%  
% % Scene two is the target plus the mask, with higher target mask
% TargetFreq = [3,4,5,6,7,8,9];
% for ii = 1:length(TargetFreq)
%     Params.freq(1) = TargetFreq(ii);
%     MaskTarget = sceneCreate('harmonic',Params);
%     MaskTarget = sceneSet(MaskTarget,'name','MaskTarget');
%     MaskTarget  = sceneAdjustLuminance(MaskTarget,meanL);
%     vcAddAndSelectObject(MaskTarget);
%     sceneWindow;
%     
%     MaskTargetXYZ = sceneGet(MaskTarget,'xyz');
%     % Go into Opponent Color space
%     % imgNOpp = imageLinearTransform(MaskTargetXYZ, colorTransformMatrix('xyz2opp', 10));
%     % And blur
%     [MaskTargetFilteredXYZ, MaskTargetFilteredOpp] = scOpponentFilter(MaskTargetXYZ,scP);
%     
%     %% Find the error xyz2 = xyz1 + E;
%     [r,c,w] = size(MaskFilteredXYZ);
%     totalE = MaskTargetFilteredXYZ - MaskFilteredXYZ;
%     method = 'ls';    % Least-squares
%     alpha = metricsMaskedError(MaskFilteredXYZ(:),totalE(:),method);
%     maskedE = alpha*MaskFilteredXYZ(:);
%     unmaskedE = totalE(:) - alpha*MaskFilteredXYZ(:);
%     
%     % Now create a new error image that differs from the original just by this
%     % decorrelated (unmasked) error
%     unmaskedXYZ2 = MaskFilteredXYZ(:) + unmaskedE(:);
%     maskedXYZ2 = MaskFilteredXYZ(:) + maskedE(:);
%     
%     % Shape it back into the image format, and display it
%     unmaskedXYZ2 = reshape(unmaskedXYZ2,r,c,w);
%     maskedXYZ2 = reshape(maskedXYZ2,r,c,w);
%     
%     %% Calculate DE
%     % What should WhiteXYZ be set to?
%     max(MaskXYZ(2,:))
%     max(MaskTargetXYZ(2,:))
%     WhiteXYZ = sceneGet(MaskTarget,'illuminantxyz'); %
%     
%      WhiteXYZ = sceneGet(MaskTarget,'illuminantxyz'); %
%     params.deltaEversion = '2000';   % Which CIELAB version
%     deltaEImage1 = scComputeDifference(MaskFilteredXYZ,MaskTargetFilteredXYZ,WhiteXYZ,params.deltaEversion);
%     deltaEImage2 = scComputeDifference(MaskFilteredXYZ,unmaskedXYZ2,WhiteXYZ,params.deltaEversion);
%     % vcNewGraphWin; imshow(deltaEImageScaled);
%     meanDE(ii) = mean(deltaEImage1(:));
%     scielab2(ii) = mean(deltaEImage2(:));
% end
% % 
% % plot(TargetFreq,meanDE,'r');
% % hold on;
% plot(TargetFreq, scielab2,'g');
% xlabel('freq');
% ylabel('DeltaE');
%     
% 
% 
% 
%    
% 
