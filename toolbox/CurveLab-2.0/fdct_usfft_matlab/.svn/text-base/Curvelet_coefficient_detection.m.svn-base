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

%%
% X = imread('bacteria.BMP');   %bacteria.BMP  calcification1.bmp  circular.PNG
% X = X(:,:,1);
% X = X(1:128,1:128);

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

% im = fdct_usfft_dispcoef(C);
% figure,imshow(im);

% % Get threshold value
% cfs =[];
% for s=1:length(C)
%   for w=1:length(C{s})
%     cfs = [cfs; abs(C{s}{w}(:))];
%   end
% end
% cfs = sort(cfs); cfs = cfs(end:-1:1);
% nb = round(pctg*length(cfs));
% cutoff = cfs(nb);
% 
% % Set small coefficients to zero
% for s=1:length(C)
%   for w=1:length(C{s})
%     C{s}{w} = C{s}{w} .* (abs(C{s}{w})>cutoff);
%   end
% end

%set the coarse resolution coefficients to zero
for a = 1:length(C{1}{1})
    for b = 1: length(C{1}{1})
        C{1}{1}(a,b) = 0;
    end
end

disp('Take inverse curvelet transform: ifdct_usfft');
tic; R_Cur = ifdct_usfft(C,0); toc;
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

%% 
% subplot(1,2,1); colormap gray; imagesc(real(X)); axis('image'); title('original image');
% subplot(1,2,2); colormap gray; imagesc(Real_Y); axis('image'); title('curvelet threshold');
subplot(2,2,1); colormap gray; imagesc(real(TestIm)); axis('image'); title('original image');
subplot(2,2,2); colormap gray; imagesc(real(R_Cur)); axis('image'); title('curvelet partial reconstruction');
subplot(2,2,3); colormap gray; imagesc(real(R_usingwavelethf)); axis('image'); title('wavelet partial reconstruction');

% %specify one curvelet
% s = 5;
% w = 1;
% [A,B] = size(C{s}{w});
% a = ceil((A+1)/2);
% b = ceil((B+1)/2);
% C{s}{w}(a,b) = 1;
% 
% %adjoint curvelet transform
% disp('Take adjoint curvelet transform: afdct_usfft');
% tic; Y = afdct_usfft(C,0); toc;
% 
% %display the curvelet
% F = ifftshift(fft2(fftshift(Y)));
% subplot(1,2,1); colormap gray; imagesc(real(Y)); axis('image'); ...
%     title('a curvelet: spatial viewpoint');
% subplot(1,2,2); colormap gray; imagesc(abs(F)); axis('image'); ...
%     title('a curvelet: frequency viewpoint');
% 
% %get parameters
% [SX,SY,FX,FY,NX,NY] = fdct_usfft_param(C);