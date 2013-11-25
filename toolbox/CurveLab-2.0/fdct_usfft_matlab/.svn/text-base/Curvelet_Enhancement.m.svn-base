clc;clear all;close all;
% I=dicomread('F:\DDSM\benign_data\case0033\C_0033_1.RIGHT_CC.dcm');
% load D:\work\svm\C_0033_1.RIGHT_CC_boundary.mat;
% I(1,1) = 0;


% figure(1),imshow(I,[]);
% hold on;plot(boundary1(2,:),boundary1(1,:));
% hold on;plot(boundary(2,:),boundary(1,:));

%% 测试图像为在原图病灶区域处切割所得512*512大小的子图像
% time = cputime;
% index2 = min(boundary(1,:));
% index1 = min(boundary(2,:))-100;
% TestIm = imcrop(I,[index1,index2,511,511]);
% % figure(2),imshow(TestIm,[]);
% TestIm = double(TestIm);
% TestIm = (TestIm-min(TestIm(:)))/(max(TestIm(:))-min(TestIm(:)));

%% Test image
X = imread('LenaCombined.jpg');   %bacteria.BMP  calcification1.bmp  circular.PNG LenaCombined.jpg
X = X(:,:,1);
% TestIm = X(1:128,1:128);
TestIm = X(288:415,238:365);

%% wavelet transform
% [Coarse,Detail] = wavedec2(TestIm,4,'db4');
% coarse_coefficients = appcoef2(Coarse,Detail,'db4',4);
% coarse_coefficients = coarse_coefficients(:);
% sizeof_coarse_coefficients = size(coarse_coefficients,1);
% 
% %set coarse resolution wavelet coefficients to zero
% Coarse(1:sizeof_coarse_coefficients) = zeros(1,sizeof_coarse_coefficients);
% 
% %Reconstraction of image
% R_usingwavelethf = waverec2(Coarse,Detail,'db4');

%% forward curvelet transform
disp('Take curvelet transform: fdct_usfft');
tic; C = fdct_usfft(TestIm,0); toc;

%% Compute norm of curvelets (Monte Carlo) 
E = cell(size(C));
 for s=1:length(C)
   E{s} = cell(size(C{s}));
   for w=1:length(C{s})
     A = C{s}{w};
%      E{s}{w} = median(abs(A(:)))/.6745;
     E{s}{w} = median(abs(A(:) - median(A(:))))/.6745; % Estimate noise level with robust estimator
   end
 end
 
%% Coefficient modified
Index_P = 0.5;Index_s = 0.6;Index_c = 3;
type = 'Mu_Cur';
Enhance_Coeff = Sigmoid_Scale_Modify(C,E,Index_P,Index_s,Index_c,type);

%% Reconstruction using modified coefficient
disp('Take inverse curvelet transform: ifdct_usfft');
tic; R_Cur = ifdct_usfft(Enhance_Coeff,0); toc;
% Real_Y = real(Y);
% % vector_ry = Real_Y(:);
% % mean_Y = mean(vector_ry);
% % std_Y = std(vector_ry);
% small_coeffs_number = find(Real_Y<=0|Real_Y>0.03);
% code = length(small_coeffs_number);
% while (code~=0)
%     Real_Y(small_coeffs_number(code)) = 0;
%     code = code-1;
% end

%% Display the result
% subplot(1,2,1); colormap gray; imagesc(real(X)); axis('image'); title('original image');
% subplot(1,2,2); colormap gray; imagesc(Real_Y); axis('image'); title('curvelet threshold');
subplot(2,2,1); colormap gray; imagesc(real(TestIm)); axis('image'); title('original image');
subplot(2,2,2); colormap gray; imagesc(real(R_Cur)); axis('image'); title('curvelet partial reconstruction');
% subplot(2,2,3); colormap gray; imagesc(real(R_usingwavelethf)); axis('image'); title('wavelet partial reconstruction');
