clear all; 
X1 = imread('circular.PNG');
X1 = X1(:,:,1);
X1 = X1(50:305,40:295);
tic; C_C = fdct_usfft(X1,0); toc;

X2 = imread('normal.PNG');
X2 = X2(:,:,1);
X2 = X2(1:256,1:256);
tic; C_N = fdct_usfft(X2,0); toc;

X3 = imread('stellate.PNG');
X3 = X3(:,:,1);
X3 = X3(1:256,1:256);
tic; C_St = fdct_usfft(X3,0); toc;