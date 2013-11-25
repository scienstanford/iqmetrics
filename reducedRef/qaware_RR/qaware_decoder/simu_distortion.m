function imwd = simu_distortion(imw, type, para, showFig)

imw = double(imw);

%% JPEG2000 compression
if (type == 1)
   imwrite(imw, (0:255)'*[1 1 1]/255, 'tmp.bmp', 'bmp');
   command = ['jp2bin\kdu_compress -i tmp.bmp -o tmp.j2k -rate ' num2str(para) ' -no_weights'];
   dos(command);
   command = ['jp2bin\kdu_expand -i tmp.j2k -o tmp_j2k.bmp'];
   dos(command);
   imwd = double(imread('tmp_j2k.bmp'));
end
%% JPEG compression
if (type == 2)
	imwrite(imw/255, 'tmp.jpg', 'quality', para);
   imwd = double(imread('tmp.jpg'));
end

%% White Gaussian noise
if (type == 3)
   sigma = para;
	imwd = imw + sigma*randn(size(imw)); 
end
%% Guassian blur
if (type == 4)
	imw = double(imw);
   g = fspecial('gaussian', [21 21], para);
   imwd = filter2(g, imw, 'same');
end
%% No distortion
if (type == 5)
   imwd = imw;
end

imwd(imwd < 0) = 0;
imwd(imwd > 255) = 255;

if (showFig)
	figure; showIm(imw, [0 255], 1);
   figure; showIm(imwd, [0 255], 1);
end

imwd = uint8(imwd);

return