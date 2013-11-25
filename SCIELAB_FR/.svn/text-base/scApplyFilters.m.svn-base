function dstImage = scApplyFilters(srcImage, filters, dimension, imgPadMethod)
% Apply the Spatial CIELAB smoothing filters to an image
%
%     scApplyFilters(oppImg, filters, dimension)
%
% The routine applies different spatial filters to the color planes of a
% source image (srcImage).
%
% The parameter filters should be a cell array with dimension equal to the
% number of color planes. If filters is not a cell, but just a single
% matrix, it will be applied to all of the color planes.
%
% The  image "srcImage" should be in RGB (Rows x Cols x Colors) format.
% Typically the srcImage is in opponent-colors space so that the first
% dimension is luminance, the other two are red-green and blue yellow.
%
% Before filtering, the function matches the sizes of the image and the 
% filters. If image size is larger than filter size, the filter will be 
% padded with zeros. However, in case that the image size is smaller than
% the filter size, we pad the image using the symmetric padding method. We 
% let the user the option to pad the image with black or white. This
% is determined by the argument 'imgPadMethod'. The user can insert '0' to
% pad with zeros or '1' to pad with ones (i.e. max value of the image).
%
%
% The parameter dimension is something I never understood, left over from a
% long time ago with XZ. We should use dimension = 2 as a default and
% ultimately get rid of it.
%
% Example:
%
%
% Copyright ImagEval Consultants, LLC, 2000.

if ieNotDefined('srcImage'),  error('src image required'); end
if ieNotDefined('filters'),   error('smoothing filters required'); end
if ieNotDefined('dimension'), dimension = 2; end
% if ieNotDefined('imgPadMethod'), imgPadMethod = 'symmetric'; end


%
if (dimension == 1), error('Why is dimension 1?'); end

[M N L] = size(srcImage);
if ~ (L == 3)
    warning('Source image is not opponent');
end

% Make sure the filters are in a cell array format
if ~iscell(filters)
    temp = filters; clear filters;
    filters=cell(1,3);
    for ii=1:L, filters{ii}=temp; end;
end;

% these varaibles are needed in order to determine the size of the padding
rowImgPad = 0;
colImgPad = 0;
rowFiltPad = 0;
colFiltPad = 0;

% Adjust the size of the filter to match the image
% This is needed so the fft2 can work properly.  Padding the filter with
% zeros should be OK because they are zero on the outside anyway. The image
% will be padded with max, zeros or symmetric padding.
fSize   = size(filters{1});
imgSize = size(srcImage(:,:,1));

% If they aren't equal, do stuff
if ~isequal(fSize,imgSize)
    
    % image size is too small - pad image
    if fSize(1) > imgSize(1) || fSize(2) > imgSize(2)
        if ieNotDefined('imgPadMethod')
            imgPadMethod = ieReadStringOrNumber('Image size smaller than filter size. This can lead to unexpected results. Please choose a padding method (symmetric/ replicate/ 0/ 1(max)):','symmetric');
        end
        % number of pixels to pad image. Max needed because only row or col
        % may be too small.
        rowImgPad = max(0,fSize(1) - imgSize(1));
        colImgPad = max(0,fSize(2) - imgSize(2));
        
        
        % If row or col is odd, we remove a row/col from the srcImage
        % This is because padding image symmetrically is preferred.
        % Otherwise the filter will shift the image data. We would rather
        % kill off a row or column in the image.
        if isodd(rowImgPad)
            srcImage = srcImage(1:(end-1),:,:);
            rowImgPad = rowImgPad + 1;
        end
        if isodd(colImgPad)
            srcImage = srcImage(:,1:(end-1),:);
            colImgPad = colImgPad + 1;
        end
    end
    % in case image is larger than filter - pad filter
    if fSize(1) < imgSize(1) || fSize(2) < imgSize(2)
        % How unequal are they?
        rowFiltPad = max(0,imgSize(1) - fSize(1));
        colFiltPad = max(0,imgSize(2) - fSize(2));

        % If row or col is odd, we remove a row/col from the srcImage
        % This is because padding filter symmetrically is preferred.
        % Otherwise the filter will shift the image data. We would rather
        % kill off a row or column in the image.
        if isodd(rowFiltPad)
            srcImage = srcImage(1:(end-1),:,:);
            rowFiltPad = rowFiltPad - 1;
        end
        if isodd(colFiltPad)
            srcImage = srcImage(:,1:(end-1),:);
            colFiltPad = colFiltPad - 1;
        end

    end
end

% Apply the filters using fft2
dstImage = zeros(size(srcImage));
for ii=1:L
    
    thisImage = srcImage(:,:,ii);
    % image needs to be padded
    if rowImgPad || colImgPad
        % imgPadMethod = 1 => pad with max image value
        if imgPadMethod == 1
            imgPadMethod = max(thisImage(:));
        end
        thisImage = padarray(thisImage,[rowImgPad/2, colImgPad/2], imgPadMethod);
    end

    thisFilter = filters{ii};
    %filter needs to be padded
    if rowFiltPad || colFiltPad
        thisFilter = padarray(thisFilter,[rowFiltPad/2, colFiltPad/2]);
    end

    % Don't we need some fft-shifting here?
    % thisImage = ones(3,3); thisImage = padarray(thisImage,[1,1]);
    % thisFilter = zeros(5,5); thisFilter(1,1) = 1;
    dstPadImg = ifftshift(real(ifft2( ...
          fft2(fftshift(thisImage)) ...
        .*fft2(fftshift(thisFilter)))));
    
    dstImage(:,:,ii) = dstPadImg(rowImgPad/2+1:end-rowImgPad/2,...
        colImgPad/2+1:end-colImgPad/2);
    % dstImage(:, :, ii) = real(ifft2(fft2(thisFilter).*fft2(thisImage)));
    % dstImage(:, :, ii) = ieConv2FFT(srcImage(:, :, ii), filters{ii}, 'same');

end;

return;
