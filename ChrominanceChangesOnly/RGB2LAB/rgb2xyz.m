function XYZ = rgb2xyz(RGB)

%RGB2XYZ converts an image from RGB format to XYZ
%
%XYZ = rgb2xyz(RGB)
%
%Output parameter: RGB must be an image having 3 channels, i.e. a 3D array.
%   Each channel of RGB should be in the range 0..255.
%Input parameter: XYZ would be an image of the same dimension.
%
%Example: To convert an image in RGB format to L*a*b* format, one can
%combine the calling to the rgb2xyz and xyz2lab functions:
%   rgbim = imread('yourImage.png');  % image in RGB format
%   labim = xyz2lab(rgb2xyz(rgbim));
%The output labim has the following ranges of values:
%  - first channel L*       0..100
%  - the second channel a*  -128..127
%  - the third channel b    -128..127
% 
%Other software packages, e.g. OpenCV, often applies a normalisation
%process so that images in the L*a*b* format can be stored as 3-channel
%greyscale images. The normalisation adopted by OpenCV is:
%  - the first channel L* is rescaled to the range 0..255. This can be done
%    easily by multiplying by 255 and dividing by 100.
%  - the second channel a* and third channel b* are added by 128 to bring
%    the range from -128..127 to 0..255.
%
%When comparing the pixel values of an L*a*b* image produced by the Matlab
%function here and those from OpenCV, ensure that the above normalisation
%procedure is taken into account.
%
%Also noted that images that are read from disk into OpenCV have their
%channels specified in reverse order in the data structure. So, to compare
%the Matlab function here with those in OpenCV, use the following OpenCV
%functions:
%   source_rgb_image = cvLoadImage('yourImage.png', 1);
%   // make sure that the 3rd parameter is CV_BGR2Lab and not CV_RGB2Lab
%   cvCvtColor(source_rgb_image, destination_lab_image, CV_BGR2Lab);
%   // the resultant L*a*b* image is created within the program (i.e. not
%   // read from disk), so the channels are in the right order:
%   cvCvtPixToPlane(destination_lab_image, l_channel, a_channel, b_channel,
%                   0);
%
%References:
%* http://www.easyrgb.com/index.php?X=MATH
%* http://en.wikipedia.org/wiki/SRGB_color_space
%
%SEE ALSO
%  xyz2rgb, xyz2lab, lab2xyz
%
%Feb 2010
%
%Copyright Du Huynh
%School of Computer Science and Software Engineering
%The University of Western Australia

RGB = RGB/255;

index1 = find(RGB > 0.04045);
index2 = find(RGB <= 0.04045);
RGB(index1) = ((RGB(index1) + 0.055) / 1.055).^2.4;
RGB(index2) = RGB(index2) / 12.92;

%Observer. = 2?бу, Illuminant = D65
XYZ = zeros(size(RGB));
M =  [
    0.4124  0.3576  0.1805;
    0.2126  0.7152  0.0722;
    0.0193  0.1192  0.9505
    ];

for i=1:3
    XYZ(:,:,i) = RGB(:,:,1)*M(i,1) + RGB(:,:,2)*M(i,2) + RGB(:,:,3)*M(i,3);
end

XYZ = XYZ*100;

end

