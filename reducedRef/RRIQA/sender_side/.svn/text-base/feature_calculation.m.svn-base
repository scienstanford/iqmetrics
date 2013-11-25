function [ggd_paras, NQ_ggd_paras] = feature_calculation(im)

level = 3;
orient = 4;
guardband = 16;
[pyr,pind] = buildSFpyr(im, level, orient-1);

wd = 2;
hd = 100;
g = fspecial('gaussian', [21 21], 3);
C = 1;

width = wd;
half_domain = hd;
domain = -half_domain:width:half_domain;
[M nbins] = size(domain);
gb = guardband;
for n = 1:orient
   band = pyrBand(pyr, pind, n+1);
   
   normfactor = filter2(g, (abs(band).^2), 'same') + C;
   nband = band./normfactor;
   nband = nband(gb+1:end-gb, gb+1:end-gb);
   [M N] = size(nband);
   lband = nband(1:M-1, 1:N-1);
   rband = nband(1:M-1, 2:N);
   dband = nband(2:M, 1:N-1);
   band1 = [lband; lband];
   band2 = [rband; dband];
   c = corrcoef(abs(band1), abs(band2));
   corr_norm = c(1, 2);
   
   band = band(gb+1:end-gb, gb+1:end-gb);
	[M N] = size(band);
   band = reshape(band, 1, M*N);
	h = hist(band, domain)./(M*N);
   ggd_paras(n, 1:2) = ggd_fit(h, width);
   
   h_hat = ggd_hist(ggd_paras(n, :), width, nbins);
   h_hat = h_hat/sum(h_hat);
   ggd_paras(n, 3) = kld(h_hat, h);
end

width = wd*4;
half_domain = hd*4;
domain = -half_domain:width:half_domain;
[M nbins] = size(domain);
gb = guardband/2;
for n = orient+1:orient*2
	band = pyrBand(pyr, pind, n+1);
   
   normfactor = filter2(g, (abs(band).^2), 'same') + C;
   nband = band./normfactor;
   nband = nband(gb+1:end-gb, gb+1:end-gb);
   [M N] = size(nband);
   lband = nband(1:M-1, 1:N-1);
   rband = nband(1:M-1, 2:N);
   dband = nband(2:M, 1:N-1);
   band1 = [lband; lband];
   band2 = [rband; dband];
   c = corrcoef(abs(band1), abs(band2));
   corr_norm = c(1, 2);
   
   band = band(gb+1:end-gb, gb+1:end-gb);
	[M N] = size(band);
	band = reshape(band, 1, M*N);
	h = hist(band, domain)./(M*N);
   ggd_paras(n, 1:2) = ggd_fit(h, width);
   
   h_hat = ggd_hist(ggd_paras(n, :), width, nbins);
   h_hat = h_hat/sum(h_hat);
   ggd_paras(n, 3) = kld(h_hat, h);
end

width = wd*16;
half_domain = hd*16;
domain = -half_domain:width:half_domain;
[M nbins] = size(domain);
gb = guardband/4;
for n = orient*2+1:orient*3
	band = pyrBand(pyr, pind, n+1);
   
   normfactor = filter2(g, (abs(band).^2), 'same') + C;
   nband = band./normfactor;
   nband = nband(gb+1:end-gb, gb+1:end-gb);
   [M N] = size(nband);
   lband = nband(1:M-1, 1:N-1);
   rband = nband(1:M-1, 2:N);
   dband = nband(2:M, 1:N-1);
   band1 = [lband; lband];
   band2 = [rband; dband];
   c = corrcoef(abs(band1), abs(band2));
   corr_norm = c(1, 2);
   
   band = band(gb+1:end-gb, gb+1:end-gb);
	[M N] = size(band);
	band = reshape(band, 1, M*N);
	h = hist(band, domain)./(M*N);
   ggd_paras(n, 1:2) = ggd_fit(h, width);
   
   h_hat = ggd_hist(ggd_paras(n, :), width, nbins);
   h_hat = h_hat/sum(h_hat);
   ggd_paras(n, 3) = kld(h_hat, h);
end

% quantization of ggd_parameters
NQ_ggd_paras = ggd_paras; % save the results without quantization

NBITS_s = 8;
NPBITS_s = 3;
NBITS_p = 8;
NBITS_d = 8;
range_p = 2^NBITS_p - 1;
range_d = 2^NBITS_d - 1;

STEP_s = 9/(2^NBITS_s);
STEP_p = 0.00625;
MIN_p = 0.05;
MAX_p = MIN_p + STEP_p*range_p;
STEP_d = 0.00025;
MIN_d = 0;
MAX_d = MIN_d + STEP_d*range_d;

t = ggd_paras(:, 1);
ep = floor(log10(t));
nt = t.*(10.^(-ep));
t = (STEP_s*round((nt-1)/STEP_s) + 1).*(10.^ep);
ggd_paras(:, 1) = t;

t = ggd_paras(:, 2);
t = round((t-MIN_p)/STEP_p)*STEP_p + MIN_p;
t(t<MIN_p) = MIN_p;
t(t>MAX_p) = MAX_p;
ggd_paras(:, 2) = t;

t = ggd_paras(:, 3);
t = round((t-MIN_d)/STEP_d)*STEP_d + MIN_d;
t(t<MIN_d) = MIN_d;
t(t>MAX_d) = MAX_d;
ggd_paras(:, 3) = t;

return