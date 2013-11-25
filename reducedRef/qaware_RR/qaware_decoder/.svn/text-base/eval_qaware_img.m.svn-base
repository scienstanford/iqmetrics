function quality = eval_qaware_img(imwd)

dwinfo = dec_wmimg(imwd);
[dinfo failure err_rate] = ecd_info(dwinfo);
if (failure)
   quality = -1;
   return
else
   dparas = feature_dequantization(dinfo);
   quality = quality_analysis(imwd, dparas);
   return
end