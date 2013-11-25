function Mout = plot_wavelet(M, Jmin, options)

% plot_wavelet - plot 2D wavelet transform stored using Mallat's ordering.
%
%   plot_wavelet(MW, Jmin, options);
%
%   'MW' is the wavelet transform (in Mallat's ordering, not lifting inplace
%       ordering).
%   'Jmin' is the minimum scale of the transform.
%
%   You can set options.style, options.edge_color, options.renormalize
%       options.line_width.
%   
%   Copyright (c) 2006 Gabriel Peyré

if nargin<2
    Jmin = 1;
end

options.null = 0;
if isfield(options, 'style')
    style = options.style;
else
    style = 'real';
end
if isfield(options, 'edge_color')
    edge_color = options.edge_color;
else
    edge_color = 'r';
end
if isfield(options, 'style')
    renormalize = options.renormalize;
else
    renormalize = 0;
end
if isfield(options, 'line_width')
    lw = options.line_width;
else
    lw = 2;
end

n = size(M,1);
Jmax = log2(n)-1;

for j=Jmin:Jmax
    qmin = double(~(j==Jmin));
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
            if q>0
                M1 = 0.3 * M1/std(M1(:));
                M1 = rescale( clamp(M1, -1,1) );
            else
                M1 = rescale(M1);
            end
        elseif strcmp(style, 'bw')
            M1 = rescale(abs(M1));
            M1 = double( M1<20/255 );
        else
            error('unknown style');    
        end
        M(selx,sely,:) = M1;
    end
    % white borders
    if 0 && edge_color>=0
        M(2^j,1:2^(j+1)) = edge_color; 
        M(1:2^(j+1),2^j) = edge_color;
    end
end

hold on;
s = [1/2 n-1/2];
imagesc(s,s,M);
axis image; axis off;
colormap gray(256);
% display axis separation
for j=Jmax:-1:Jmin;
    x = [0 2^(j+1)];
    y = [2^j 2^j];
    h = plot(x,y, edge_color);
    set(h, 'LineWidth', lw);
    y = [0 2^(j+1)];
    x = [2^j 2^j];
    h = plot(x,y, 'r');
    set(h, 'LineWidth', lw);
end
% plot overall square
x = [0 0 n n 0];
y = [0 n n 0 0];
h = plot(x,y,edge_color);
set(h, 'LineWidth', lw*1.2);
hold off;
axis ij;

if nargout>1
    Mout = M;
end