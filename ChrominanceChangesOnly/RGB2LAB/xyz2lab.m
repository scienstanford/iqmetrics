function LAB = xyz2lab(XYZ)

%XYZ2LAB converts an image from XYZ format to CIE-L*a*b*
%
%LAB = xyz2LAB(XYZ)
%
%Input parameter: XYZ must be an image having 3 channels, i.e. a 3D array.
%Output parameter: LAB would be an image of the same dimension.
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
%  lab2xyz, xyz2rgb, rgb2xyz
%
%Feb 2010
%
%Copyright Du Huynh
%School of Computer Science and Software Engineering
%The University of Western Australia

%ref_X =  95.047   Observer= 2?бу, Illuminant= D65
%ref_Y = 100.000
%ref_Z = 108.883

ref_XYZ = [95.047; 100.000; 108.883];

for i=1:3
    XYZ(:,:,i) = XYZ(:,:,i) / ref_XYZ(i);   
end

index1 = find(XYZ > 0.008856);
index2 = find(XYZ <= 0.008856);
XYZ(index1) = XYZ(index1).^(1/3);
XYZ(index2) = (7.787 * XYZ(index2) ) + ( 16 / 116 );

% CIE-L*a*b*
LAB = zeros(size(XYZ));
LAB(:,:,1) = ( 116 * XYZ(:,:,2) ) - 16;
LAB(:,:,2) = 500 * ( XYZ(:,:,1) - XYZ(:,:,2) );
LAB(:,:,3) = 200 * ( XYZ(:,:,2) - XYZ(:,:,3) );
end

