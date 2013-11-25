% clear all;
% info = dicominfo('E:\MATLAB6p5\work\DICOMIMAGE\face.dcm');
% X = dicomread(info);
% 
% % clear all;
% % load woman;
% 
% 
% figure;
% imshow(X, []);
% title('Original');
% pixval;
% 
% X=double(X);
% YuantuZuizhi =max(max(X)');
% XReverse=YuantuZuizhi-X;
% figure;
% imshow(XReverse,[]);
% title('XReverse');
% pixval;
% 
% XReverse=uint16(XReverse);
% XReverseHisteq=Histeq(XReverse,64);
% figure;
% imshow(XReverseHisteq,[]);
% title('XReverseHisteq');
% pixval;
% 
% figure;
% imshow(XReverse,[]);
% title('XReversequzheng');
% pixval;
% 
% imhist(XReverse,64);title('XReverse');
% figure; imhist(XReverseHisteq,64);title('XReverseHisteq');
function Lvbo=wavelateenhancement(X)

% clear all;
% load woman;
% figure;
% imshow(X,[]);
% title('Original');
% 
% YuantuZuizhi =max(max(X)');


% I = rgb2gray(I);
% 
% 
% [ca,ch,cv,cd] = dwt2( X,'db1','mode','sym');
% [ca,ch,cv,cd] = dwt2(X,'db4');
%  
% 
% figure;%SUBPLOT(2,2,1);
% imshow(ca, []);
% title('Decompose Approximation Coefficients');
% 
% figure;% SUBPLOT(2,2,2);
% imshow(ch, []);
% title('Decompose Horizontal Coefficients');
% 
% figure;% SUBPLOT(2,2,3);
% imshow(cv, []);
% title('Decompose Vertical Coefficients');
% 
% figure;% SUBPLOT(2,2,4);
% imshow(cd, []);
% title('Decompose Diagonal Coefficients');

X = imread('bacteria.BMP');
% X = X(:,:,1);
X = X(1:128,1:128);
figure,imshow(X);
% 分解
n = 4;
[C,S]= wavedec2(X,n,'db4');%C :分解的近似和细节信息A（n），H（n），V（n），D（n），H（n－1），V（n－1），D（n－1）...
                           %S :各层图像的大小信息S（1，：），S（i，：）...
cAn = appcoef2(C,S,'db4',n);
%[H,V,D]=detcoef2('all',C,S,n);
cAn1 = cAn(:);% row_wise 排列，形成一列
m= size(cAn1,1);
w = C(m+1:length(C));%取出后面的细节分量

%系数调整
max0 = max(abs(w));%细节图像的最值
% max0 =YuantuZuizhi;
g = 3;%门限内映射斜率
% t = 1;
t = 0.25;
%调整门限计算
d = max0*t; 
% d1=YuantuZuizhi*0.25;

g0 = (g-1)*d;%门限外调整值

i = find(w>d);
w(i) = w(i)+g0;
% w(i) = w(i);

j = find(w<-d);
w(j) = w(j)-g0;
% w(i) = w(i);

k = find(abs(w)<=d);
w(k)=w(k)*g;    

%重构
C(m+1:length(C))=w;
Chonggou = waverec2(C,S,'db4');
Shuchu = uint8(Chonggou);
figure;imshow(Shuchu,[]);
Lvbo = wiener2(Shuchu,[3 3]);
figure;imshow(Lvbo,[]);
% title('Output:g=3 t=0.25');

HisteqTu = histeq(Lvbo);
figure;imshow(HisteqTu,[]);
% title('Output:g=3 t=0.025,HisteqLater');
  