function info = feature_quantization(paras)

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

selected_bands = 1:size(paras,1);
len = length(selected_bands);

% encode features into binary array
s = paras(selected_bands, 1)';%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ep = floor(log10(s));
nt = s.*(10.^(-ep));
qs = abs(nt - 1)/STEP_s;
qsep(1:3)  =  ep(1:3) - (-6);
qsep(4:7)  =  ep(4:7) - (-5);
qsep(8:11) =  ep(8:11) - (-4);
% qsep(11:14)=  ep(11:14) - (-2);


   
p = paras(selected_bands, 2)';
qp = round((p - MIN_p)/STEP_p);
   
d = paras(selected_bands, 3)';
qd = round((d - MIN_d)/STEP_d);
   
len_s = NBITS_s*len;
len_sp = NPBITS_s*len;
len_p = NBITS_p*len;
len_d =  NBITS_d*len;
bins = reshape(dec2bin(qs, NBITS_s)', 1, len_s);
binsp = reshape(dec2bin(qsep, NPBITS_s)', 1, len_sp);
binp = reshape(dec2bin(qp, NBITS_p)', 1, len_p);
bind = reshape(dec2bin(qd, NBITS_d)', 1, len_d);
infostr = [bins binsp binp bind];
   
% information bits

for j=1:length(infostr)
   
   info(j) = str2num(infostr(j));
end

return