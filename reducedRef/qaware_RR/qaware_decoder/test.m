clear all;
close all;

imw = imread('qa_img.tif');

%imwd = simu_distortion(imw, 1, 0.5, 1);
%distortion = eval_qaware_img(imwd)

%imwd = simu_distortion(imw, 2, 18, 1);
%distortion = eval_qaware_img(imwd)

imwd = simu_distortion(imw, 3, 15, 1);
distortion = eval_qaware_img(imwd)

%imwd = simu_distortion(imw, 4, 1.5, 1);
%distortion = eval_qaware_img(imwd)
