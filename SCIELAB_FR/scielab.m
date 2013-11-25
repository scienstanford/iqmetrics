function [result, params, xyz1, xyz2] = scielab(image1,image2,whitePt,params)
%Spatial CIELAB (S-CIELAB) metric
%
%  [deltaEImage, params, xyz1, xyz2] = scielab(image1,image2,whitePt,[params])
%
% The Spatial-CIELAB difference image (error map) compares image1
% and image2.  The metric is compatible with CIELAB over constant portions
% of the image.  It adds spatial blurring consistent with measurements of
% human space-color sensitivity to account for the image spatial structure.
%
% image1 and image2: 3-D images in XYZ or LMS format. 
% whitePt:     a cell  array containing the white points of the two images. 
% params:      a structure containing several variables used in the
%              calculation. The entires are updated and can be returned 
%              by the routine. 
%
%  params.
%       sampPerDeg = How many samples per degree of visual angle in the image.
%                    If the image is, say, 5 deg, and contains 128 samples,
%                    then this parameter is 512/2. 
%                    The default is 224 for historical reasons.  In
%                    general, the code should be improved to work well at
%                    low sample rates.
%       filters    = filters used in spatial blurring of opponent channels.
%                    If these are present, then the filters are used.
%                    Otherwise, new filters are created. They will be
%                    returned in params to save time in the next call.
%       filterSize = usually equal to sampPerDeg. 
%       imageFormat= Data format of the input image.
%             'xyz2', 'xyz10', 'lms2', 'lms10';
%       deltaEversion  = which version of CIELAB.  (Default is 2000)
%              Earlier options permit '1976', '1994' and '2000'. 
%              Added for special ISET analyses, we allow a request for
%              CIELAB 2000 'chroma','hue', or 'luminance' component errors.
%              These are always calculated using CIELAB 2000.
%
% The routine is divided into three main sub-routines that perform the
% key computations.
%
%   scPrepareFilter     -- Prepare the spatial blurring filters
%   scOpponentFilter    -- Apply the filters in opponent color space
%   scComputeDifference -- Compute the resulting CIELAB differences
%
%
% Example:
%   See scielabExample.m for a full description
%
%   scielab(image1,image2,whitePt,params);  % deltaE 2000 CIELAB version
%
%   params.deltaEversion = '1976';
%   scielab(image1,image2,whitePt,params);  % deltaE 1976 CIELAB version
%
%   params.deltaEversion = 'hue';           % Just the hue error
%   scielab(image1,image2,whitePt,params);
%   params.deltaEversion = 'chrominance';   % Just the chrominance error
%   scielab(image1,image2,whitePt,params);
%
%   params.deltaEversion = '2000';   % Which CIELAB version
%   params.sampPerDeg = sampPerDeg;  % Sets up the viewing distance
%   params.imageFormat = imageformat; %
%   params.filterSize = sampPerDeg;
%   params.filters = [];             % Not precomputed
%  [errorImage,params] = scielab(img1LMS, img2LMS, whitePt, params);
% 
% Copyright ImagEval Consultants, LLC, 2003.

% Initialize parameters
defaultSampPerDeg = 224;

if ieNotDefined('image1'), errordlg('Scielab requires image1'); end
if ieNotDefined('image2'), errordlg('Scielab requires image2'); end
if ieNotDefined('whitePt'), errordlg('Scielab requires a white point or white point cell array'); end
if ieNotDefined('params')
    params.deltaEversion = '2000';
    params.sampPerDeg = defaultSampPerDeg;   % A big number
    params.filterSize = defaultSampPerDeg;   % 1 deg
    params.filters = [];
    params.imageFormat = 'xyz10';
    % Check if the input images are 1-D or 2-D.  I don't understand this.
    if (size(image1,1)>1 && size(image1,2)>3),  params.dimension = 2;
    else                                        params.dimension = 1;
    end
end
   
% The white point used in the CIELAB calculation.  This is expected to be a
% cell array containing a white point for each image.  If it is just a
% vector, then we convert it to a cell array.
if ~iscell(whitePt)
    tmp{1} = whitePt; tmp{2} = whitePt;
    whitePt = tmp;
end

% If the image and data format are not XYZ, we must change the white points
% to be in XYZ 
if strncmp(params.imageFormat,'lms',3)
    for i=1:2
        % Why are we changing to opp here?
        whitePt{i} = changeColorSpace(whitePt{i}, cmatrix('lms2opp'));
        if strncmp(params.imageFormat,'xyz10',5) || ...
                strncmp(params.imageFormat,'lms10',5),
            xyztype = 10;
        else xyztype = 2;
        end
        % Why is this opp2xyz and not lms2xyz from the beginning?  Is there
        % no cmatrix for lms2xyz?
        whitePt{i} = changeColorSpace(whitePt{i}, cmatrix('opp2xyz', xyztype));
    end;
end

% These are the filters for spatial blurring.  They can take a
% while to create (and we should speed that up).   
if isempty(params.filters)
    [params.filters, params.support] = scPrepareFilters(params); 
end
% figure; imagesc(params.filters{2})
% Filter the image in opponent-colors space starting from lms or xyz.  The
% returned image is in XYZ.
xyz1 = scOpponentFilter(image1,params);
% figure; imagesc(xyz1(:,:,2))
xyz2 = scOpponentFilter(image2,params);

result = scComputeDifference(xyz1,xyz2,whitePt,params.deltaEversion);

% modified by WEN
% [m,n] = size(result);
% D0=0.8;
% k =sum(sum(abs(result)))/m*n;
% % load k_contourlet;
% % k_contourlet = [k_contourlet k];
% % save k_contourlet k_contourlet;
% % k =sum(1-((NC2-NC1)./NC1).^2);
% % quality = D1/(D1+log2(1+k/D0));
% quality=1/(1+log2(1+k'/D0));

return;

