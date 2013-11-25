function imw = enc_wmimg(im, einfo, showFig)
% embed bitstream "einfo" into "im"

% load embedding key and coding parameters

load key_pattern
load code_paras

% information bits
info1 = einfo(1:Nbits1);
info2 = einfo(Nbits1+1:Nbits1+Nbits2);
info3 = einfo(Nbits1+Nbits2+1:Nbits1+Nbits2+Nbits3);

% wavelet transform
im = double(im);
[pyr, pind] = buildWpyr(im, level);

%%%%% band 1
Nbits = Nbits1;
info = info1;
Key = Key1(1:Nbits);
wbandind = wbandind1;
band = pyrBand(pyr, pind, wbandind);
[M N] = size(band);
qband = band(2:M-1, 2:N-1);
qband = reshape(qband, prod(size(qband)), 1);
qband(Key) = round(qband(Key)/delta)*delta;
qband(Key) = qband(Key) + (info*2 - 1)*delta/4;
qband = reshape(qband, M-2, N-2);
band(2:M-1, 2:N-1) = qband;
pyr(pyrBandIndices(pind, wbandind)) = reshape(band, prod(size(band)), 1);

%%%%% band 2
Nbits = Nbits2;
info = info2;
Key = Key2(1:Nbits);
wbandind = wbandind2;
band = pyrBand(pyr, pind, wbandind);
[M N] = size(band);
qband = band(2:M-1, 2:N-1);
qband = reshape(qband, prod(size(qband)), 1);
qband(Key) = round(qband(Key)/delta)*delta;
qband(Key) = qband(Key) + (info*2 - 1)*delta/4;
qband = reshape(qband, M-2, N-2);
band(2:M-1, 2:N-1) = qband;
pyr(pyrBandIndices(pind, wbandind)) = reshape(band, prod(size(band)), 1);

%%%%% band 3
Nbits = Nbits3;
info = info3;
Key = Key3(1:Nbits);
wbandind = wbandind3;
band = pyrBand(pyr, pind, wbandind);
[M N] = size(band);
qband = band(2:M-1, 2:N-1);
qband = reshape(qband, prod(size(qband)), 1);
qband(Key) = round(qband(Key)/delta)*delta;
qband(Key) = qband(Key) + (info*2 - 1)*delta/4;
qband = reshape(qband, M-2, N-2);
band(2:M-1, 2:N-1) = qband;
pyr(pyrBandIndices(pind, wbandind)) = reshape(band, prod(size(band)), 1);

% reconstruct image
imw = reconWpyr(pyr, pind);
imw(imw<0) = 0;
imw(imw>255) = 255;

if (showFig)
	figure; showIm(im, [0 255], 1);
   figure; showIm(imw, [0 255], 1);
end

imw = uint8(imw);

return