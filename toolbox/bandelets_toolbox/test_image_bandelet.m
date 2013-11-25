clear;close all;clc;
load lena;
options.use_single_qt = 0;
%  options.use_single_qt = 1;
 M = size(lena,1);
 N = size(lena,2);

 Jmin = 1;%%level of dec is 3
 j_min = 2;
 j_max = 4;
 T = 15;
 s = Inf;
 
 [QT,Theta] = compute_wavelet_quadtree(lena,Jmin,T,j_min,j_max,s, options);
 
 %%Bandelet Tranform
 [MB,Rbandg] = perform_wavelet_bandelet_transform(lena,Jmin,QT,Theta,1, options);
