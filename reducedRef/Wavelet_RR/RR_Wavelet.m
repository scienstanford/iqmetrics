function distortion=RR_Wavelet(im_reference,im_distorted,alpha)
%得到基于WBCT变换的部分参考的评价数据
%alpha为阈值权值
path('F:\code\RR_Wavelet\sender_side',path);
path('F:\code\RR_Wavelet\receiver_side',path);
path('F:\code\Wavelet9_7_toolbox',path);
NC1= sender_feature_extraction(im_reference,alpha);%%确定参考图像里面大于阈值的系数的比例
distortion = receiver_distortion_measure(im_distorted, NC1);%得到经过降质图像的质量测度 Q值
rmpath('F:\code\Wavelet9_7_toolbox');
rmpath( 'F:\code\RR_Wavelet\receiver_side');
rmpath('F:\code\RR_Wavelet\sender_side');