function Q = DCTune_FR(img_ref, img_dis)

imwrite(img_ref,'ref.ppm');
imwrite(img_dis,'dis.ppm');

cmd = 'F:\code\FR_DCTune\dctune2.0 -error ref.ppm dis.ppm';

[status, result] = system(cmd);

delete('ref.ppm');
delete('dis.ppm');

if status
    beg = findstr(result, 'Perceptual Error:') + 17;
    q = result(beg:end);
else
    error('Failed to execute the DCTune command');
end

Q = str2num(q);
    