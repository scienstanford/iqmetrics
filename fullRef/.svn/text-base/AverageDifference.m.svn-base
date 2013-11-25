%% Program for Average Difference Calculation



function AD = AverageDifference(origImg, distImg)

origImg = double(origImg);
distImg = double(distImg);

[M N] = size(origImg);
error = origImg - distImg;

AD = sum(sum(error)) / (M * N);