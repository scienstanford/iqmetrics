function d = kld(hist1, hist2)

hist1 = hist1/sum(hist1);
hist2 = hist2/sum(hist2);
ind = hist1>0 & hist2>0;

d = sum(hist1(ind).*log(hist1(ind)./hist2(ind)))/sum(hist1(ind));
d = abs(d);

return