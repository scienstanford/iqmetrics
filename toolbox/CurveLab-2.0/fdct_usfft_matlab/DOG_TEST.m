clc;close all;clear all;

% I=dicomread('E:\cancer_data\case1530\A_1530_1.RIGHT_MLO.dcm');
% load A_1530_1.RIGHT_MLO.1_boundary.mat;
% boundary1=boundary;
% load A_1530_1.RIGHT_MLO.2_boundary.mat;

% I=dicomread('F:\DDSM\cancer_data\case3395\B_3395_1.RIGHT_CC.dcm');
% load B_3395_1.RIGHT_CC_boundary.mat;
% boundary1=boundary;
% load B_3515_1.LEFT_CC.2_boundary.mat;

% I=dicomread('E:\cancer_data\case4125\D_4125_1.LEFT_MLO.dcm');
% load D_4125_1.LEFT_MLO_boundary.mat;

% I=dicomread('H:\benign_data\case1442\A_1442_1.LEFT_MLO.dcm');
% load A_1442_1.LEFT_MLO_boundary.mat;
% boundary1=boundary;
% load C_0248_1.LEFT_MLO.2_boundary.mat;

% I=dicomread('H:\benign_data\case0033\C_0033_1.RIGHT_MLO.dcm');
% load C_0033_1.RIGHT_MLO_boundary.mat;

% I=dicomread('H:\benign_data\case3100\B_3100_1.LEFT_MLO.dcm');
% load B_3100_1.LEFT_MLO_boundary.mat;

% I=double(I);
% I=mat2gray(I);
% I(1,1) = 0;
% figure,imshow(I,[]);
% hold on;plot(boundary1(2,:),boundary1(1,:));
% hold on;plot(boundary(2,:),boundary(1,:));

load CropImage_3395.mat;
[r,c]=size(I2);
E_Image=I2(1:4:r,1:4:c);
figure,imshow(E_Image);
h = fspecial('gaussian',[3,3],20);