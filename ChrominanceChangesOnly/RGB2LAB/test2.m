%  RGB2XYZ
clear all
addpath(genpath('D:\tools\matlab\workspace\iset-4.449\iset-4.0'));
RGB = imread('E:\2012\Error Map\TID\I03.bmp');
 RGB1 = double(imread('E:\2012\Error Map\TID\I03.bmp'));

XYZ1 = rgb2xyz(RGB1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X1 = XYZ1(:,:,1); 
X3= Arnold(X1,1,0);
[m,n] = size(X1);
% x1 = X1-20;  %  0.1~1
x1 = X1+randn([m,n])*5; % the different x can influence the result of SSIM
Y1 = XYZ1(:,:,2); 
Y3 = Y1/100;
%y1 = Y1+randn([m,n])*200;
%y1 = Y1-60;
Z1 = XYZ1(:,:,3);
% z1 = Z1+randn([m,n])*100; % <120
 z1 = Z1+10;
XYZ2(:,:,1)=X3;
XYZ2(:,:,2)=Y1; 
XYZ2(:,:,3)=Z1;

plot(X1,X3,'k-');



% %%%%%%%%%%%%%%%%%%%%%%%%%%
% CIE Lab space
% 
% lab1 = xyz2lab(XYZ1);
% L1 = lab1(:,:,1);
% [m,n] = size(L1);
% %l1 = L1+randn([m,n])*200;
% %l1 = L1+20; 
% a1 = lab1(:,:,2);
% 
% A1 = a1+randn([m,n])*200;
%   %A1 = a1+20; 
% b1 = lab1(:,:,3);
% % B1 = b1+randn([m,n])*20;  % <200
%  B1 = b1+60;     % < 32
% 
% lab2(:,:,1)=L1;
% lab2(:,:,2)=a1; 
% lab2(:,:,3)=B1;
% 
% XYZ2= lab2xyz(lab2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 RGB2 =xyz2rgb(XYZ2);

XYZ3 = rgb2xyz(RGB2);

X2 = XYZ3(:,:,1); 
Y2 = XYZ3(:,:,2);
Z2 = XYZ3(:,:,3); 

[PSNRY,error_mapY] = PeakSignaltoNoiseRatio(Y1, Y2);
[mssimY, ssim_mapY] = SSIM_FR(Y1, Y2);
[PSNRX,error_mapX] = PeakSignaltoNoiseRatio(X1, X2);
[mssimX, ssim_mapX] = SSIM_FR(X1, X2);
[PSNRZ,error_mapZ] = PeakSignaltoNoiseRatio(Z1, Z2);
[mssimZ, ssim_mapZ] = SSIM_FR(Z1, Z2);

 % plot(X1,X2,Y1,Y2,Z1,Z2);
% figure(1);
% plot(Y1,Y2,'r-');
% xlabel('Y1');
% ylabel('Y2');
% figure(2);
% plot(X1,X2,'c-');
% xlabel('X1');
% ylabel('X2');
% figure(3);
% plot(Z1,Z2,'k-');
% xlabel('Z1');
% ylabel('Z2');
% hold on;


% figure;
% plot(Y1,Y2,'-k');
% xlabel('Y1');
% ylabel('Y2');
% figure;
% plot(X1,X2,'ok');
% xlabel('X1');
% ylabel('X2');
% figure;
% plot(Z1,Z2,'--k');
% xlabel('Z1');
% ylabel('Z2');



RGB5=uint8(RGB2);
figure, imshow(RGB5);
 title('');
% imwrite(RGB2,'C:\Users\mr\Desktop\color distortion\XYZ\y_20.bmp');
 RGB3 = RGB2GRAY(RGB1/255);
 RGB4 = RGB2GRAY(RGB2/255);
%  figure, imshow(RGB3);
%   figure, imshow(RGB4);
% [mssim1, ssim_map] = SSIM_FR(RGB3, L1);
 [mssim, ssim_map] = SSIM_FR(RGB3, RGB4);
 [PSNR,error_map] = PeakSignaltoNoiseRatio(RGB3, RGB4);
%  plot(RGB3, RGB4,'k-');
%  xlabel('');
%  ylabel('RGB3/RGB4')
 plot(Y3,RGB3,'k-');
 xlabel('');
 ylabel('Y1/RGB')