%% Program for Structural Content Calculation


function SC = StructuralContent(origImg, distImg)

origImg = double(origImg);
distImg = double(distImg);

SC = sum(sum(origImg .* origImg)) / sum(sum(distImg .* distImg));