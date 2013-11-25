function quality = quality_analysis(im, NC1)

NC2 = rr_feature_calculation(im,NC1);
% k= kld(NC1,NC2);
D0=0.8;
k =sum(abs(NC1(2:end) - NC2(2:end)));
% load k_contourlet;
% k_contourlet = [k_contourlet k];
% save k_contourlet k_contourlet;
% k =sum(1-((NC2-NC1)./NC1).^2);
% quality = D1/(D1+log2(1+k/D0));
quality=1/(1+log2(1+k/D0));
% quality=1/(1+1000*log10(1+k/D0));
% quality=1/(1+log(1+k/D0)/log(1.5));
% quality=k;
return
