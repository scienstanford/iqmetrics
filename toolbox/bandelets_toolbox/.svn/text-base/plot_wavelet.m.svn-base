function plot_wavelet(M, L, style, edge_color, renormalize)

% plot_wavelet - plot 2D wavelet transform stored using Mallat's ordering.
%
%   plot_wavelet(WT, L, style, edge_color, renormalize);
%
%   'WT' is the wavelet transform (in Mallat's ordering, not lifting inplace
%       ordering).
%   'Jmin' is the minimum scale of the transform.
%   
%   Copyright (c) 2003 Gabriel Peyré

if nargin<2
    L = 1;
end

if nargin<3 || isempty(style)
    style = 'real';
end

if nargin<4
    edge_color = 0;
end

if nargin<5
    renormalize = 0;
end

n = length(M);

for j=L:log2(n)-1
    qmin = ~(j==L);
    for q=qmin:3
        [selx,sely] = compute_quadrant_selection(j,q);
        M1 = M(selx,sely,:);
        
        if strcmp(style, 'abs')
            M1 = rescale( abs(M1) );

            if renormalize
                M1 = 1.5*M1;
                I = find(abs(M1)>1);
                M1(I) = M1(I)./abs(M1(I));
            end

        elseif strcmp(style, 'absinv')
            
            M1 = rescale(abs(M1));
            if renormalize
                M1 = 2*M1;
                I = find(abs(M1)>1);
                M1(I) = M1(I)./abs(M1(I));
            end

            M1 = rescale( -abs(M1) );
        elseif strcmp(style, 'real')
            M1 = rescale(M1);
        elseif strcmp(style, 'bw')
            M1 = rescale(abs(M1));
            M1 = double( M1<20/255 );
        else
            error('unknown style');    
        end
        
        
        M(selx,sely,:) = M1;
    end
    % white borders
    if edge_color>=0
        M(2^j,1:2^(j+1)) = edge_color; 
        M(1:2^(j+1),2^j) = edge_color;
    end
end

imagesc(M);
axis image;
axis off;