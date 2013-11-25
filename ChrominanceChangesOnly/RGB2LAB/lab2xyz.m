function XYZ = lab2xyz(LAB)

%LAB2XYZ converts an image from CIE-L*a*b* format to XYZ
%
%XYZ = lab2xyz(LAB)
%
%Input parameter: L*a*b* must be an image having 3 channels, i.e. a 3D array.
%   The first (L) channel must have values in the range 0..100; the second
%   and third (a* and b*) channels must have values in the range -128 to
%   127.
%Output parameter: XYZ would be an image of the same dimension.
%
%The XYZ format is often used as an intermediate format for conversion from
%RGB to L*a*b* and vice versa.
%
%Example: To convert an image in RGB format to L*a*b* format, one can
%combine the calling to the rgb2xyz and xyz2lab functions:
%   rgbim = imread('yourImage.png');  % image in RGB format
%   labim = xyz2lab(rgb2xyz(rgbim));
%The output labim has the following ranges of values:
%  - first channel L*       0..100
%  - the second channel a*  -128..127
%  - the third channel b*   -128..127
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
%  xyz2lab, xyz2rgb, rgb2xyz
%
%Feb 2010
%
%Copyright Du Huynh
%School of Computer Science and Software Engineering
%The University of Western Australia


XYZ = zeros(size(LAB));

XYZ(:,:,2) = ( LAB(:,:,1) + 16 ) / 116;
XYZ(:,:,1) = LAB(:,:,2) / 500 + XYZ(:,:,2);
XYZ(:,:,3) = XYZ(:,:,2) - LAB(:,:,3) / 200;

index1 = find(XYZ.^3 > 0.008856);
index2 = find(XYZ.^3 <= 0.008856);

XYZ(index1) = XYZ(index1).^3;
XYZ(index2) = ( XYZ(index2) - 16 / 116 ) / 7.787;

ref_XYZ = [95.047; 100.000; 108.883];
for i=1:3
    XYZ(:,:,i) = XYZ(:,:,i) * ref_XYZ(i);
end
end

%Observer= 2?бу, Illuminant= D65