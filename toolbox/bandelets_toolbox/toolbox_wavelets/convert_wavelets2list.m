function y = convert_wavelets2list(x, Jmin, Jmax, pack_imagettes);

% convert_wavelets2list - convert a wavelet transform into a cell array.
%
%	y = convert_wavelets2list(x, Jmin, Jmax, pack_imagettes);
%
%   If x is a wavelet transform, then x is a cell array (a list)
%       of cell array. y{1} is a cell array containing
%       H/V/D quadrant of the transform at finest scale (Jmax), 
%       y{2} contains the quadrant at scale Jmax-1, etc.
%       y{end} contains the low scale decomposition.
%   If x is a call array, then it compute the inverse decomposition.
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<4
    pack_imagettes = 1;
end

if ~iscell(x)
    d = size(x,3);
    if nargin<3
        Jmax = log2(size(x,1))-1;
    end
    if nargin<2
        error('You must provide Jmin');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % turn the wavelets coefficient into a cell array
    y = {};
    for j=Jmax:-1:Jmin
        % concatenate the 3 imagettes
        M_list = {};
        for q=1:3
            for m=1:d
                [selx,sely] = compute_quadrant_selection(j,q);
                if pack_imagettes
                    M_list{q+(m-1)*d} = x(selx,sely,m);
                else
                    M_list{m} = x(selx,sely,m);
                end
            end
            if ~pack_imagettes
                y{end+1} = M_list;
                M_list = {};
            end
        end
        if pack_imagettes
            y{end+1} = M_list;
        end
    end
    % add the low scale image
    M_list = {};
    for m=1:d
        M_list{m} = x(1:2^Jmin,1:2^Jmin,m);
    end
    y{end+1} = M_list;
else
    Jmax = log2(size(x{1}{1},1));
    if pack_imagettes
        d = length(x{1})/3;
        Jmin = Jmax - length(x)+2;
    else
        d = length(x{1});
        Jmin = Jmax - (length(x)-1)/3+1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % turn the cell array into a list
    n = 2^(Jmax+1);
    y = zeros(n,n,d);
    k=0;
    for j=Jmax:-1:Jmin
        if pack_imagettes
            k = k+1;
            M_list = x{k};
        end
        % concatenate the 3 imagettes
        for q=1:3
            if ~pack_imagettes
                k = k+1;
                M_list = x{k};
            end
            for m=1:d
                [selx,sely] = compute_quadrant_selection(j,q);
                if pack_imagettes
                    y(selx,sely,m) = M_list{q+(m-1)*d};
                else
                    y(selx,sely,m) = M_list{m};
                end
            end
        end
    end
    % add the low scale image
    M_list = x{end};
    for m=1:d
        y(1:2^Jmin,1:2^Jmin,m) = M_list{m};
    end
end