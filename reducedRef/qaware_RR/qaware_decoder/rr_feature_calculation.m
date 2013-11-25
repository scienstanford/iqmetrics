function klds = rr_feature_calculation(im, paras)

level = 3;
orient = 4;
guardband = 16;
[pyr,pind] = buildSFpyr(double(im), level, orient-1);

wd = 2;
hd = 100;

klds = zeros(1, level*orient);
ind = [1 3 6 8 9 11];

width = wd;
half_domain = hd;
domain = -half_domain:width:half_domain;
[M nbins] = size(domain);
gb = guardband;
for n = ind(ind<=orient)
	band = pyrBand(pyr, pind, n+1);
   band = band(gb+1:end-gb, gb+1:end-gb);
	[M N] = size(band);
	band = reshape(band, 1, M*N);
	h = hist(band, domain)./(M*N);
   h_ref = ggd_hist(paras(n, :), width, nbins);
   h_ref = h_ref/sum(h_ref);
	klds(n) = kld(h_ref, h);
end

width = wd*4;
half_domain = hd*4;
domain = -half_domain:width:half_domain;
[M nbins] = size(domain);
gb = guardband/2;
for n = ind(ind>orient & ind<=(2*orient))
	band = pyrBand(pyr, pind, n+1);
   band = band(gb+1:end-gb, gb+1:end-gb);
	[M N] = size(band);
	band = reshape(band, 1, M*N);
	h = hist(band, domain)./(M*N);
   h_ref = ggd_hist(paras(n, :), width, nbins);
   h_ref = h_ref/sum(h_ref);
	klds(n) = kld(h_ref, h);
end

width = wd*16;
half_domain = hd*16;
domain = -half_domain:width:half_domain;
[M nbins] = size(domain);
gb = guardband/4;
for n = ind(ind>(2*orient) & ind<=(3*orient))
	band = pyrBand(pyr, pind, n+1);
   band = band(gb+1:end-gb, gb+1:end-gb);
	[M N] = size(band);
	band = reshape(band, 1, M*N);
	h = hist(band, domain)./(M*N);
   h_ref = ggd_hist(paras(n, :), width, nbins);
   h_ref = h_ref/sum(h_ref);
	klds(n) = kld(h_ref, h);
end

return