function paras = feature_dequantization(info)

NBITS_s = 8;
NPBITS_s = 3;
NBITS_p = 8;
NBITS_d = 8;
range_p = 2^NBITS_p - 1;
range_d = 2^NBITS_d - 1;

STEP_s = 9/(2^NBITS_s);
STEP_p = 0.00625;
MIN_p = 0.05;
MAX_p = MIN_p + STEP_p*range_p;
STEP_d = 0.00025;
MIN_d = 0;
MAX_d = MIN_d + STEP_d*range_d;

selected_bands = 1:11;

len = length(selected_bands);

% recover features from binary array
for j=1:length(info)
   rinfo(j) = num2str(info(j));
end
len_s = NBITS_s*len;
len_sp = NPBITS_s*len;
len_p = NBITS_p*len;
len_d =  NBITS_d*len;
rbins = rinfo(1:len_s);
rbinsp = rinfo(len_s+1:len_s+len_sp);
rbinp = rinfo(len_s+len_sp+1: len_s+len_sp+len_p);
rbind = rinfo(len_s+len_sp+len_p+1:end);
rqs = bin2dec(reshape(rbins, NBITS_s, len)')';%二进制变成十进制
rqsep = bin2dec(reshape(rbinsp, NPBITS_s, len)')';
rqp = bin2dec(reshape(rbinp, NBITS_p, len)')';
rqd = bin2dec(reshape(rbind, NBITS_d, len)')';
rs(1:3) = (rqs(1:3)*STEP_s + 1).*10.^(rqsep(1:3)-6);
rs(4:7) = (rqs(4:7)*STEP_s + 1).*10.^(rqsep(4:7)-5);
rs(8:11) = (rqs(8:11)*STEP_s + 1).*10.^(rqsep(8:11)-4);
% rs(11:14) = (rqs(11:14)*STEP_s + 1).*10.^(rqsep(11:14)-2);
rp = rqp*STEP_p + MIN_p;
rd = rqd*STEP_d + MIN_d;

% write back to paras
paras = zeros(11, 3);
paras(selected_bands, 1) = rs';
paras(selected_bands, 2) = rp';
paras(selected_bands, 3) = rd';

return