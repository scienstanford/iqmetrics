%得到经过降质图像的质量测度  Q值
function quality = quality_analysis(im, NC1)

NC2 = rr_feature_calculation(im,NC1);%确定降质图像里各个方向子带中幅值大于一定阈值系数在该子带内所占的比例
% k= kld(NC1,NC2);
D0=0.1;
k =sum(abs(NC1(2:end) - NC2(2:end)));
% load k_WBCT;
% k_WBCT = [k_WBCT k];
% save k_WBCT k_WBCT;
% k =sum(1-((NC2-NC1)./NC1).^2);
% quality = D1/(D1+log2(1+k/D0));
% quality=1/(1+log2(1+k/D0));
quality=1/(1+log10(1+k/D0));
% quality=1/(1+log(1+k/D0)/log(1.5));
% quality=k;
return
