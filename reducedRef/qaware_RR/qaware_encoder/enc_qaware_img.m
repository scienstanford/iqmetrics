function imw = enc_qaware_img(im)

paras = feature_calculation(im);
info = feature_quantization(paras);
einfo = ecc_info(info);
imw = enc_wmimg(im, einfo, 1);

return