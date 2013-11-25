function bins = ggd_hist(paras, width, nbins)
	samp_dens = 30;

	half_range = width*nbins/2;
	step = width/samp_dens;
	x = -half_range+step/2:step:half_range;
   p = ggd(x, paras(1), paras(2));
   [M N] = size(x);
	p = reshape(p, samp_dens, N/samp_dens);
   bins = sum(p);
   if (bins==0)
      bins = ones(1, nbins);
   end
return
