function quality = quality_analysis(im, rparas)

klds = rr_feature_calculation(im, rparas);

ind = [1 3 6 8 9 11];
D0 = 0.001;
kld = sum(abs(klds(ind) - rparas(ind, 3)'));
quality = log2(1 + kld/D0);

return
