function dinfo = dec_wmimg(imwd)

% load coding parameters and embedding key
load key_pattern
load code_paras

imwd = double(imwd);
[pyr, pind] = buildWpyr(imwd, level);

% decode information bits

%%%%% band1
Nbits = Nbits1;
Key = Key1(1:Nbits);
wbandind = wbandind1;
band = pyrBand(pyr, pind, wbandind);
[M N] = size(band);
qband = band(2:M-1, 2:N-1);
qband = reshape(qband, prod(size(qband)), 1);
dband(Key) = 1 - mod(floor(qband(Key)/(delta/2)), 2);
dinfo1 = dband(Key);

%%%%% band2
Nbits = Nbits2;
Key = Key2(1:Nbits);
wbandind = wbandind2;
band = pyrBand(pyr, pind, wbandind);
[M N] = size(band);
qband = band(2:M-1, 2:N-1);
qband = reshape(qband, prod(size(qband)), 1);
dband(Key) = 1 - mod(floor(qband(Key)/(delta/2)), 2);
dinfo2 = dband(Key);

%%%%% band3
Nbits = Nbits3;
Key = Key3(1:Nbits);
wbandind = wbandind3;
band = pyrBand(pyr, pind, wbandind);
[M N] = size(band);
qband = band(2:M-1, 2:N-1);
qband = reshape(qband, prod(size(qband)), 1);
dband(Key) = 1 - mod(floor(qband(Key)/(delta/2)), 2);
dinfo3 = dband(Key);

% save decoded bits
dinfo = [dinfo1'; dinfo2'; dinfo3'];

return