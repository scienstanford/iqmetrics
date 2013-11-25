function distortion=RR_Contourlet(im_reference,im_distorted,alpha)
path('F:\code\RR_Contourlet\sender_side',path);
path('F:\code\RR_Contourlet\receiver_side',path);
path('F:\code\contourlet_toolbox',path);
NC1= sender_feature_extraction(im_reference,alpha);
distortion = receiver_distortion_measure(im_distorted, NC1);
rmpath('F:\code\contourlet_toolbox');
rmpath( 'F:\code\RR_Contourlet\receiver_side');
rmpath('F:\code\RR_Contourlet\sender_side');