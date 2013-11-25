function m = evaluate_nbr_bits(M,T)

% evaluate_nbr_bits - evaluate the number of bits to code an image
%
%   m = evaluate_nbr_bits(M)
%
%   nbr_bit is a lower bound (Shanon entropy bound) on the number 
%       of bits needed to code the wavelet image.
%
%   Works also with cell array of matrix.
%
%   Copyright (c) 2004 Gabriel Peyré

if iscell(M)
    m = 0;
    for i=1:length(M)
        m = m + evaluate_nbr_bits(M{i},T);
    end
    return;
end
    
[y, nbr_bits, nbr_bits_tot, M] = perform_quantization(M, T);

M = M(:);
n = length(M);
h = []; % histogram
for p = unique(M)';
    h = [h, length(find(M==p))];
end
h = h/n;
m = - n * sum( h.*log2(h) );
