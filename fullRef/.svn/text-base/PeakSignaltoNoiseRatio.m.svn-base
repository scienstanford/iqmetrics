%% Program for Peak Signal to Noise Ratio Calculation



function PSNR = PSNR(origImg, distImg)

origImg = double(origImg);
distImg = double(distImg);

[M N] = size(origImg);
error = origImg - distImg;
MSE = sum(sum(error .* error)) / (M * N);

if(MSE > 0)
    PSNR = 10*log(255*255/MSE) / log(10);
else
    PSNR = 99;
end