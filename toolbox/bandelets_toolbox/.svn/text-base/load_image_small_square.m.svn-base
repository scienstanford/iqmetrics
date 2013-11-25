function M = load_image_small_square(name, p, options)

% load_image_small_square - extract a sub-square from an image.
%
%   M = load_image_small_square(name, p, options);
%
%   p is the size of the small square.
%   name is the name of the undelying image, please see 'help load_image'
%   options can override default features.
%
%   Copyright (c) 2005 Gabriel Peyré

options.null = 0;

if isfield(options, 'use_wavelets')
    use_wavelets = options.use_wavelets;
else
    use_wavelets = 1;
end

if isfield(options, 'use_filtering')
    use_filtering = options.use_filtering;
else
    use_filtering = 0;
end

if isfield(options, 'wavelet_scale')
    wavelet_scale = options.wavelet_scale;
else
    wavelet_scale = 0;
end

if isfield(options, 'n')
    n = options.n;
else
    n = 256;
end

% set up the ROI for extraction
c = [0.5,0.5];        % default roi center
switch lower(name)
    case 'lena'
        c = [0.8022, 0.6767];
    case 'lena_bis'
        c = [0.875,0.665];
    case 'quarterdisk'
        c = [0.38,0.72];
    case 'real_wall_1024'
        c = [0.51,0.57];
    case 'peppers'
        c = [0.821,0.306];
    case 'peppers_line'
        c = [0.821,0.306];     % line
    case 'peppers_curve'
        c = [0.21,0.136];      % curve
    case 'barb'
        c = [0.478,0.62];
    case 'barb_left_pant'
        c = [0.662,0.388];     % left pant
    case 'real_light'
        c = [0.59,0.32];        
    case 'line'
        options.theta = 1/sqrt(2);
        options.eta = 0.5-options.theta*0.5;
end

% load the original image
M = load_image(name, n, options);

% perform a 2D wavelet transform
if use_wavelets
    Jmin = 1;
    Jmax = log2(n)-1;
    j = Jmax - wavelet_scale;
    q = 1;
    % use 7-9 wavelet transform
    options.wavelet_vm = 4;
    MW = perform_wavelet_transform(M,Jmin,1,options);
    % select sub-square
    [selx,sely] = compute_quadrant_selection(j,q); 
    MW = MW(selx,sely);
elseif use_filtering>0
    sigma = use_filtering;
    MW = conv2(M, ones(sigma)/sigma^2, 'same');    
else
    MW = M;
end

% crop sub-image
n = length(MW);
cn = floor(c*(n-1))+1;
p1 = floor(p/2); p2 = ceil(p/2);
M = MW( cn(1)-p1+1:cn(1)+p2, cn(2)-p1+1:cn(2)+p2 );