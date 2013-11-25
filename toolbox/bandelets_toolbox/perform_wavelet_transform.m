function MW = perform_wavelet_transform(M,Jmin,dir, options)

% perform_wavelet_transform - a wavelet transform
%
%   MW = perform_wavelet_transform(M,Jmin,dir,options);
%
%   You can replace this code with your favorite
%   wavelet transform.
%
%   MW must be stored in "Mallat's order", 
%   ie ie has the same size as M, with various
%   scales and orientation packed from up-left 
%   to bottom right.
%
%   Copyright (c) 2005 Gabriel Peyr?

options.null = 0;
if nargin<2
    Jmin = 1;
end
if nargin<3
    dir=1;
end

if size(M,3)>1
    MW = M;
    for i=1:size(M,3)
        MW(:,:,i) = perform_wavelet_transform(M(:,:,i),Jmin, dir, options);
    end
    return;
end


% number of dimension
ndim = length(size(M));
if ndim==2 && ( size(M,2)==1 || size(M,1)==1 )
    ndim=1;
end

global wavelet_vm;
if isfield(options, 'wavelet_vm')
    wavelet_vm = options.wavelet_vm;
else
    if isempty(wavelet_vm) || ndim==2
        wavelet_vm = 4;
    end
end

% test if wavelab is available
if ~exist('MakeBSFilter')
    wavelet_vm = 0; % use haar transform, sorrrrry
end

% retrieve the 7-9 CDF biorthogonal filters
if wavelet_vm>1
    [qmf,dqmf] = MakeBSFilter( 'CDF', [wavelet_vm,wavelet_vm] );
end

% compute biorthogonal wavelet transform
if ndim==1
    if wavelet_vm>1
        if dir==1
            MW = FWT_SBS(M,Jmin,qmf,dqmf);
        else
            MW = IWT_SBS(M,Jmin,qmf,dqmf);
        end
    else
        MW = perform_haar_transform_slow(M,Jmin,dir);
    end
else
    if exist('qmf') && wavelet_vm>0
        if dir==1
            MW = FWT2_SBS(M,Jmin,qmf,dqmf);
        else
            MW = IWT2_SBS(M,Jmin,qmf,dqmf);
        end
    else
        % no wavelab: use mex haar transform instead (it uses Jmin==1)
        if dir==-1
            % makes conversion Mallat's order -> inplace order
            M = reorder_coefs(M,0,-1);
        end
        MW = perform_haar_transform(M,Jmin,dir);
        if dir==1
            % makes conversion inplace order -> Mallat's order
            MW = reorder_coefs(MW,0,1);
        end
    end
end