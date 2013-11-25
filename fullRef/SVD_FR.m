function [scaMeasure, graMeasure] = SVD_IQA(origImg, distImg)
%% SVD-based Image Quality Measure
% This method outputs a graphical & numerical image quality measure
% based on Singular Value Decomposition.

%% For details on the implementation, please refer
% Aleksandr Shnayderman, Alexander Gusev, and Ahmet M. Eskicioglu,
% "An SVD-Based Grayscale Image Quality Measure for Local and Global Assessment",
% IEEE TRANSACTIONS ON IMAGE PROCESSING, VOL. 15, NO. 2, FEBRUARY 2006.
%% Parameters
% img1          -   Input Reference Gray Image
% img2          -   Input Distorted Gray Image
% blkSize       -   Window size for block processing
% graMeasure    -   Graphical Image quality measure
% scaMeasure    -   Numerical Image quality measure
%% Sample
% img1=imread('hats.bmp');
% figure;subplot(121);
% imshow(img1);title('Orginal Image');
% 
% img2= imread('hats_JPEG.bmp');
% subplot(122);
% imshow(img2);title('Distorted Image');
% 
% blkSize = 8;
% 
% [scaMeasure,graMeasure] = SVD_IQA(img1, img2, blkSize);
% 
% figure;imshow(graMeasure);title('Graphical Measure');
% disp('Numerical Measure');
% disp(scaMeasure);
%% Initization

blkSize = 8;
k = size(origImg, 1);
% ///begain change
if nargin<3
    blkSize=11;
end
% ///end change
blkx = blkSize;
blky = blkSize;

blockwise1 = MatDec(origImg,blkx);
blockwise2 = MatDec(distImg,blkx);
[blkx blky imgx imgy] = size(blockwise1);
graMeasure = zeros(imgx,imgy);
blockwise1 = double(blockwise1);
blockwise2 = double(blockwise2);

%% Run this function

for i=1:imgx
    for j=1:imgy
        temp_in_image = blockwise1(:,:,i,j);temp_in_image=temp_in_image(:);
        original_img = reshape(temp_in_image,blkx,blky);
        temp_dist_image = blockwise2(:,:,i,j);temp_dist_image=temp_dist_image(:);
        distorted_img = reshape(temp_dist_image,blkx,blky);
        graMeasure(i,j) = sqrt(sum((svd(original_img)-svd(distorted_img)).^2));
    end
end

graMeasure = round((graMeasure/max(max(graMeasure)))*255);
graMeasure = uint8(graMeasure);

scaMeasure = sum(sum(abs(graMeasure-median(median(graMeasure)))))/((k/blkx).^2);



%% This function decomposes an image into blocks.
%% Parameters
% inImg         -   Input Gray Image
% blkSize       -   Window size for block processing
% out           -   Output 4 dimensional matrix with blocks.
%
function out = MatDec(inImg, blkSize)

[m,n]=size(inImg);

r3=m/blkSize;
c3=n/blkSize;
q4=0;
q1=0;

for i=1:r3
    for j=1:c3
        for s=1:blkSize
            for k=1:blkSize
                p3=s+q4;
                q2=k+q1;
                out(s,k,i,j)=inImg(p3,q2);
            end
        end
        q1=q1+blkSize;
    end
    q4=q4+blkSize;q1=0;
end

