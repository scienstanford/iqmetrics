function [y, nbr_bits, nbr_bits_tot, y_quant] = perform_quantization(x, T)

% perform_quantization - perform a nearly uniform quantization of the signal.
%
%   [y, nbr_bits, nbr_bits_tot, y_quant] = perform_quantization(x, T);
%
%   The quantizer is defined by y=Q_T(x) where:
%       Q_T(x) = 0    if  |x|<T
%       Q_T(x) = sign(x) * ([x/T]+0.5)*T      where [.]=floor
%   (i.e. a nearly uniform quantizer with twice larger zero bin).
%
%   y is the quantified value
%   nbr_bit as same size as x and store the number of bits needed to code
%       each entry of x
%   nbr_bits_tot is just sum(nbr_bits_tot(:))
%   y_quant is the signed token representing each entry of y.
%
%   Copyright (c) 2004 Gabriel Peyré

if iscell(x)
    nbr_bits_tot = 0;
    for i=1:length(x)
        [y{i}, nbr_bits{i}, nbr_bits_t] = perform_quantization(x{i}, T);
        nbr_bits_tot = nbr_bits_tot + nbr_bits_t;
    end
    return;
end

I = find(abs(x)<T);
nbr_bits = floor(abs(x)/T);

y_quant = sign(x).*nbr_bits;

y = sign(x).*(nbr_bits+0.5)*T;
y(I) = 0;

nbr_bits = nbr_bits+1;
nbr_bits(I) = 0;

nbr_bits_tot = sum(nbr_bits(:));