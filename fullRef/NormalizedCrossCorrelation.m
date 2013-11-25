%% Program for Normalized Cross Correlation Calculation


function NK = NormalizedCrossCorrelation(origImg, distImg)

origImg = double(origImg);
distImg = double(distImg);

NK = sum(sum(origImg .* distImg)) / sum(sum(origImg .* origImg));