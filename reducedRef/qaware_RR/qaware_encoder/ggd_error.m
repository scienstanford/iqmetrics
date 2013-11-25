function error = ggd_error(paras, hist, width)

[M nbins] = size(hist);
h = ggd_hist(paras, width, nbins);
h = h/sum(h);
%error = mean((h - hist).^2);
error = kld(h, hist);

return