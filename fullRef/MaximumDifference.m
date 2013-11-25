%% Program for Maximum Difference Calculation



function MD = MaximumDifference(origImg, distImg)

origImg = double(origImg);
distImg = double(distImg);

error = origImg - distImg;

MD = max(max(error));