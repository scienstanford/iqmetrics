function features = sender_feature_extraction(im)

[paras, nq_paras] = feature_calculation(im);
features = feature_quantization(paras);

return