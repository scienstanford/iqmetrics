%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A code for MUSICA algorithm which developed by AGFA, who get the patent
% in June 30, 1992.The modyfication of coeffcients in this algorithm came 
% from patent and Jing Xv's paper.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear all;close all;

I=dicomread('F:\DDSM\benign_data\case0033\C_0033_1.RIGHT_CC.dcm');
load D:\work\svm\C_0033_1.RIGHT_CC_boundary.mat;

I(1,1) = 0;
% figure(1),imshow(I,[]);
% hold on;plot(boundary1(2,:),boundary1(1,:));
% hold on;plot(boundary(2,:),boundary(1,:));

%% 测试图像为在原图病灶区域处切割所得512*512大小的子图像
time = cputime;
index2 = min(boundary(1,:));
index1 = min(boundary(2,:))-100;
TestIm = imcrop(I,[index1,index2,511,511]);
% figure(2),imshow(TestIm,[]);
TestIm = double(TestIm);
TestIm = (TestIm-min(TestIm(:)))/(max(TestIm(:))-min(TestIm(:)));
time = cputime-time

%% 确定分解层数
scale_index = log(size(TestIm,1))/log(2)-4;

%% Laplacian decomposition using Gaussian filters and scale_index levels
pfilt = 'GF_LP';
LP_Coeff_EL = MultiL_LPLayer(TestIm, pfilt, scale_index);

%% Display output
% figure(1)
% colormap(gray)
% nr = floor(sqrt(scale_index+1));
% nc = ceil((scale_index+1)/nr);
% for l = 1:scale_index+1
%     subplot(nr, nc, l); 
%     imageshow(LP_Coeff_EL{l});
% end
clear I boundary;
%% Perfect Reconstruction
% LP_RC_Im = MultiR_LPLayer(LP_Coeff_EL, pfilt);

%% Enhancement from MUSICA patent
% type = 'Mu_P';
% Index_P = 0.7;
% Enhance_Coeff = Sigmoid_Scale_Modify(LP_Coeff_EL,Index_P,type);
% Enhance_MUSICA_Im = MultiR_LPLayer(Enhance_Coeff, pfilt);

%% Enhancement from Mingqing Li(MUSICA I)
Index_P = 0.7;
type = 'Mu_L';
mini_bound = 3;maxi_bound = 5;
ample_modcoe = 1.7;ample_dyacoe = 1.4;
Enhance_Coeff = Sigmoid_Scale_Modify(LP_Coeff_EL,Index_P,type);
Enhance_MUSICA_Im = MultiE_LPLayer(Enhance_Coeff, pfilt,mini_bound,maxi_bound,ample_modcoe,ample_dyacoe);
% Enhance_MUSICA_Im = MultiR_LPLayer(Enhance_Coeff, pfilt);

%% Show perfect reconstruction
figure(2)
colormap(gray)
subplot(1,2,1), imshow(TestIm, []);
subplot(1,2,2), imshow(Enhance_MUSICA_Im, []);
% subplot(1,2,1), imageshow(TestIm, [0, 1]);
% subplot(1,2,2), imageshow(Enhance_MUSICA_Im, [0, 1]);
% title(sprintf('SNR = %.2f dB', SNR(x, xr)))
%%
% [histgram,maxV,minV] = Hist_of_Image(I);
% figure(3),plot(minV:maxV,histgram(:));
% [step,totalleft,totalright,CropImage,CropBinaryImage]=extragraysplit(I);

