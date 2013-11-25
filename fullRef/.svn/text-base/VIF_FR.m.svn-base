function vif= VIF_FR(origImg,distImg)


% -----------COPYRIGHT NOTICE STARTS WITH THIS LINE------------
% Copyright (c) 2005 The University of Texas at Austin
% All rights reserved.
% 
% Permission is hereby granted, without written agreement and without license or royalty fees, to use, copy, 
% modify, and distribute this code (the source files) and its documentation for
% any purpose, provided that the copyright notice in its entirety appear in all copies of this code, and the 
% original source of this code, Laboratory for Image and Video Engineering (LIVE, http://live.ece.utexas.edu)
% at the University of Texas at Austin (UT Austin, 
% http://www.utexas.edu), is acknowledged in any publication that reports research using this code. The research
% is to be cited in the bibliography as:
% 
% H. R. Sheikh and A. C. Bovik, "Image Information and Visual Quality", IEEE Transactions on 
% Image Processing, (to appear).
% 
% IN NO EVENT SHALL THE UNIVERSITY OF TEXAS AT AUSTIN BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, 
% OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OF THIS DATABASE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF TEXAS
% AT AUSTIN HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
% THE UNIVERSITY OF TEXAS AT AUSTIN SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE DATABASE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS,
% AND THE UNIVERSITY OF TEXAS AT AUSTIN HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
% 
% -----------COPYRIGHT NOTICE ENDS WITH THIS LINE------------
%
%This is an implementation of the algorithm for calculating the
%Visual Information Fidelity (VIF) measure (may also be known as the Sheikh
%-Bovik Index) between two images. Please refer
%to the following paper:
%
%H. R. Sheikh and A. C. Bovik "Image Information and Visual Quality"
%IEEE Transactios on Image Processing, in publication, May 2005.
%Download manuscript draft from http://live.ece.utexas.edu in the
%Publications link
%
%This implementation is slightly differnet from the one used to report
%results in the paper above. The modification have to do with using more
%subands than those used in the paper, better handling of image boundaries,
%and a window that automatically resizes itself based on the scale.
%
%Report bugfixes and comments to hamid.sheikh@ieee.org
%
%----------------------------------------------------------------------
% Prerequisites: The Steerable Pyramid toolbox. Available at
% http://www.cns.nyu.edu/~lcv/software.html
%
%Input : (1) img1: The reference image
%        (2) img2: The distorted image (order is important)
%
%Output: (1) VIF te visual information fidelity measure between the two images

%Default Usage:
%   Given 2 test images img1 and img2, whose dynamic range is 0-255
%
%   vif = vifvec(img1, img2);
%
%Advanced Usage:
%   Users may want to modify the parameters in the code. 
%   (1) Modify sigma_nsq to find tune for your image dataset.
%   (2) MxM is the block size that denotes the size of a vector used in the
%   GSM model.
%   (3) subbands included in the computation
%========================================================================


path(path,'D:\Tools\Matlab\workspace\Image Quality Metrics\toolbox\matlabPyrTools');



M=3;
subbands=[4 7 10 13 16 19 22 25];
sigma_nsq=0.4;

% Do wavelet decomposition. This requires the Steerable Pyramid. You can
% use your own wavelet as long as the cell arrays org and dist contain
% corresponding subbands from the reference and the distorted images
% respectively.
[pyr,pind] = buildSpyr(origImg, 4, 'sp5Filters', 'reflect1'); % compute transform
org=ind2wtree(pyr,pind); % convert to cell array
[pyr,pind] = buildSpyr(distImg, 4, 'sp5Filters', 'reflect1');
dist=ind2wtree(pyr,pind);

% calculate the parameters of the distortion channel
[g_all,vv_all]=vifsub_est_m(org,dist,subbands,M);

% calculate the parameters of the reference image
[ssarr, larr, cuarr]=refparams_vecgsm(org,subbands,M);

% reorder subbands. This is needed since the outputs of the above functions
% are not in the same order
vvtemp=cell(1,max(subbands));
ggtemp=vvtemp;
for(kk=1:length(subbands))
    vvtemp{subbands(kk)}=vv_all{kk};
    ggtemp{subbands(kk)}=g_all{kk};
end


% compute reference and distorted image information from each subband
for i=1:length(subbands)
    sub=subbands(i);
    g=ggtemp{sub};
    vv=vvtemp{sub};
    ss=ssarr{sub};
    lambda = larr(sub,:);, 
    cu=cuarr{sub};

    % how many eigenvalues to sum over. default is all.
    neigvals=length(lambda);
    
    % compute the size of the window used in the distortion channel estimation, and use it to calculate the offset from subband borders
    % we do this to avoid all coefficients that may suffer from boundary
    % effects
    lev=ceil((sub-1)/6);
    winsize=2^lev+1; offset=(winsize-1)/2;
    offset=ceil(offset/M);
    
    
    % select only valid portion of the output.
    g=g(offset+1:end-offset,offset+1:end-offset);
    vv=vv(offset+1:end-offset,offset+1:end-offset);
    ss=ss(offset+1:end-offset,offset+1:end-offset);
    
    
    %VIF
    temp1=0; temp2=0;
    for j=1:length(lambda)
        temp1=temp1+sum(sum((log2(1+g.*g.*ss.*lambda(j)./(vv+sigma_nsq))))); % distorted image information for the i'th subband
        temp2=temp2+sum(sum((log2(1+ss.*lambda(j)./(sigma_nsq))))); % reference image information
    end
    num(i)=temp1;
    den(i)=temp2;
    
end

% compuate VIF
vif=sum(num)./sum(den);

%%
function [g_all, vv_all]= vifsub_est_m(org,dist, subbands, M)

% uses convolution for determining the parameters of the distortion channel
% Called by vifvec.m

tol = 1e-15; % tolernace for zero variance. Variance below this is set to zero, and zero is set to this value to avoid numerical issues.


for i=1:length(subbands)
    sub=subbands(i);
    y=org{sub};
    yn=dist{sub};

    % compute the size of the window used in the distortion channel estimation
    lev=ceil((sub-1)/6);
    winsize=2^lev+1; offset=(winsize-1)/2;
    win = ones(winsize);
    
    % force subband size to be multiple of M
    newsize=floor(size(y)./M)*M;
    y=y(1:newsize(1),1:newsize(2));
    yn=yn(1:newsize(1),1:newsize(2));

    % Correlation with downsampling. This is faster than downsampling after
    % computing full correlation.
    winstep=[M M];
    winstart=[1 1].*floor(M/2)+1;
    winstop=size(y)-ceil(M/2)+1;
    
    % mean
    mean_x = corrDn(y,win/sum(win(:)),'reflect1',winstep, winstart,winstop);
    mean_y = corrDn(yn,win/sum(win(:)),'reflect1',winstep, winstart,winstop);
    % cov
    cov_xy = corrDn(y.*yn, win, 'reflect1',winstep, winstart,winstop) - sum(win(:)).*mean_x.*mean_y;
    % var
    ss_x = corrDn(y.^2,win, 'reflect1',winstep, winstart,winstop) - sum(win(:)).*mean_x.^2;
    ss_y = corrDn(yn.^2,win, 'reflect1',winstep, winstart,winstop) - sum(win(:)).*mean_y.^2;

    
    % get rid of numerical problems, very small negative numbers, or very
    % small positive numbers, or other theoretical impossibilities.
    ss_x(ss_x<0)=0;
    ss_y(ss_y<0)=0;
   
    % Regression 
    g = cov_xy./(ss_x+tol);
    
    % Variance of error in regression
    vv = (ss_y - g.*cov_xy)/(sum(win(:)));
    
    % get rid of numerical problems, very small negative numbers, or very
    % small positive numbers, or other theoretical impossibilities.
    g (ss_x < tol) = 0;
    vv (ss_x < tol) = ss_y (ss_x < tol);
    ss_x(ss_x<tol)=0;
    
    g (ss_y < tol) = 0;
    vv (ss_y < tol) = 0;
    
    % constrain g to be non-negative. 
    vv(g<0)=ss_y(g<0);
    g(g<0)=0;
    
    % take care of numerical errors, vv could be very small negative
    vv( vv <= tol) = tol;
    
    g_all{i}=g;
    vv_all{i}=vv;
    
end

%%
function [ssarr, l_arr, cu_arr]=refparams_vecgsm(org,subands,M)

%This function computes the parameters of the reference image. This is
%called by vifvec.m.

for i=1:length(subands);
    sub=subands(i);
    y=org{sub};
    
    sizey=floor(size(y)./M)*M; % crop  to exact multiple size
    y=y(1:sizey(1),1:sizey(2));
    
    
    % Collect MxM blocks. Rearrange each block into an
    % M^2 dimensional vector and collect all such vectors.
    % Collece ALL possible MXM blocks (even those overlapping) from the subband
    temp=[];
    for j=1:M
        for k=1:M
            temp=cat(1,temp,reshape(y(k:end-(M-k), j:end-(M-j)),1,[]));
        end
    end
    
    % estimate mean and covariance
    mcu=mean(temp')';
    cu=((temp-repmat(mcu,1,size(temp,2)))*(temp-repmat(mcu,1,size(temp,2)))')./size(temp,2); % covariance matrix for U
    
    % Collect MxM blocks as above. Use ONLY non-overlapping blocks to
    % calculate the S field
    temp=[];
    for j=1:M
        for k=1:M
            temp=cat(1,temp,reshape(y(k:M:end, j:M:end),1,[]));
        end
    end

    % Calculate the S field
    ss=(inv(cu)*temp);
    ss=sum(ss.*temp)./(M*M);
    ss=reshape(ss,sizey/M);
    
    % Eigen-decomposition
    [v,d]=eig(cu);
    l_arr(sub,:)=diag(d)';
    
    % rearrange for output
    ssarr{sub}=ss;
    temp=0;
    d=diag(d);
    cu_arr{sub}=cu;
end    

%%
function  wtree = ind2wtree(pyr, ind)

%this function is called by vifvec.m
% converts the output of Eero Simoncelli's pyramid routines into subbands in a cell array
C=pyr;
S=ind;

offset=0;
numsubs=size(ind,1);
for i=1:numsubs
    wtree{numsubs-i+1}=reshape(C(offset+1:offset+prod(S(i,:))), S(i,1),S(i,2));
    offset=offset+prod(S(i,:));
end

