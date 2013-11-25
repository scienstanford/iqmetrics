function [iwssim iwmse iwpsnr]= IW_SSIM_MSE_PSNR_FR(origImg,distImg, iw_flag, Nsc, K, L, weight, win, blk_size, parent, sigma_nsq)

MAX_PSNR = 1000;

if (nargin < 2 | nargin > 11)
    iwssim = -Inf;
    iwmse = -Inf;
    iwpsnr = -Inf;
    disp(['Input argument error']);
    return;
end

if (size(origImg) ~= size(distImg))
    iwssim = -Inf;
    iwmse = -Inf;
    iwpsnr = -Inf;
    disp(['Input argument error: images should have the same size']);
    return;
end

if (~exist('iw_flag'))
   iw_flag = 1;
end

if (~exist('Nsc'))
   Nsc = 5;
end

if (~exist('K'))
   K = [0.01 0.03]; % default from [Wang et al, IEEE-TIP 2004]
end

if (~exist('L'))
   L = 255;
end

if (~exist('weight'))
   weight = [0.0448 0.2856 0.3001 0.2363 0.1333]; % default from [Wang et al, IEEE-Asilomar 2003]
end

if (exist('win'))
    winSz = size(win);
    if (winSz(1) ~= winSz(2))
        iwssim = -Inf;
        iwmse = -Inf;
        iwpsnr = -Inf;
        disp(['Window size error']);
        return;
    end
    win_size = winSz(1);
else
    win_size = 11; % default from [Wang et al, IEEE-TIP 2004]
    gwin_std = 1.5; % default from [Wang et al, IEEE-TIP 2004]
    win = fspecial('gaussian', win_size, gwin_std); % default from [Wang et al, IEEE-TIP 2004]
    win = win/sum(win(:));
end

if (~exist('blk_size'))
   blk_size = 3; % spatial neighborhood size
end

if (~exist('parent'))
   parent = 1; % include parent neighbor
end

if (~exist('sigma_nsq'))
   sigma_nsq = 0.4; % default from [Sheikh & Bovik, IEEE-TIP 2006]
end

blSzX = blk_size; 
blSzY = blk_size;
Nsc = min(5,Nsc);
weight = weight(1:Nsc);
weight = weight/sum(weight);
if (min(size(origImg)) < win_size*(2^(Nsc-1)))
    iwssim = -Inf;
    iwmse = -Inf;
    iwpsnr = -Inf;
    disp(['Image size too small for ' num2str(Nsc) ' scale IW-SSIM evaluation.']);
    return;
end
bound = ceil((win_size-1)/2);
bound1 = bound - floor((blSzX-1)/2);

[pyro,pind]= buildLpyr(double(origImg),Nsc);
[pyrd,pind]= buildLpyr(double(distImg),Nsc);

[cs_map l_map se_map] = scale_quality_maps(pyro,pyrd,pind,Nsc,K,L,win);
if (iw_flag)
    iw_map = info_content_weight_map(pyro,pyrd,pind,Nsc,parent,blSzX,blSzY,sigma_nsq);
end

for s = 1:Nsc
    se = se_map{s};
    cs = cs_map{s};
    if (s==Nsc)
        cs = cs.*l_map;
    end
    if (iw_flag)
        if (s<Nsc)
            iw = iw_map{s};
            iw = iw(bound1+1:end-bound1, bound1+1:end-bound1);
        else
            iw = ones(size(cs));
        end
        se = se(bound+1:end-bound, bound+1:end-bound);        
        wmcs(s) = sum(sum(cs.*iw))/sum(sum(iw));
        wmse(s) = sum(sum(se.*iw))/sum(sum(iw));
    else
        wmcs(s) = mean2(cs);
        wmse(s) = mean2(se);
    end
end

iwmse = prod(wmse(1:Nsc).^weight(1:Nsc));
iwpsnr = min(10*log10(255^2/iwmse), MAX_PSNR);
iwssim = prod(wmcs(1:Nsc).^weight(1:Nsc));

%% contrast-structure similarity map and squared error map for each scale, and luminance similarity map for the coarsest scale
function [cs_map l_map se_map]= scale_quality_maps(pyro,pyrd,pind,Nsc,K,L,win)
    if (nargin < 3 | nargin > 7)
        cs_map = -Inf;
        l_map = -Inf;
    	se_map = -Inf;
        disp(['Input argument error']);
        return;
    end
    if (~exist('Nsc'))
        Nsc = size(pind,1);
    end
    if (~exist('K'))
        K = [0.01 0.03]; % default from [Wang et al, IEEE-TIP 2004]
    end
    if (~exist('L'))
        L = 255;
    end
    if (~exist('win'))
        win = fspecial('gaussian', 11, 1.5); % default from [Wang et al, IEEE-TIP 2004]
        win = win/sum(win(:));
    end

    pyro = real(pyro);
    pyrd = real(pyrd);
    C1 = (K(1)*L)^2;
    C2 = (K(2)*L)^2;

    for i=1:Nsc
        y = pyrBand(pyro, pind, i);
        yn = pyrBand(pyrd, pind, i);
        se_map{i} = (y-yn).^2;
        
        mu1 = filter2(win, y, 'valid');
        mu2 = filter2(win, yn, 'valid');
        sigma12 = filter2(win, y.*yn, 'valid') - mu1.*mu2;
        sigma1_sq  = filter2(win, y.^2, 'valid') - mu1.^2;
        sigma2_sq  = filter2(win, yn.^2, 'valid') - mu2.^2;
        sigma1_sq = max(0, sigma1_sq);
        sigma2_sq = max(0, sigma2_sq);
        cs_map{i} = (2*sigma12 + C2)./(sigma1_sq + sigma2_sq + C2);
        
        if (i == Nsc)
            l_map = (2*mu1.*mu2 + C1)./(mu1.^2 + mu2.^2 + C1);
        end
    end
    
%% compute information content weight map for Scale 1 to Nsc-1
function [iw_map]= info_content_weight_map(pyro,pyrd,pind,Nsc,parent,blSzX,blSzY,sigma_nsq)
tol = 1e-15;

if (~exist('Nsc'))
   Nsc = size(pind, 1);
end
if (~exist('parent'))
   parent = 1; % include parent neighbor
end
if (~exist('blSzX'))
   blSzX = 3;
end
if (~exist('blSzY'))
   blSzY = 3;
end
if (~exist('sigma_nsq'))
   sigma_nsq = 0.4; % default from [Sheikh & Bovik, IEEE-TIP 2006]
end
Nband = Nsc-1;    
win = ones(3); % window for estimating gain factor g
win = win/sum(sum(win));

pyro = real(pyro);
pyrd = real(pyrd);
for nband=1:Nband
    y = pyrBand(pyro, pind, nband);
    yn = pyrBand(pyrd, pind, nband);
    
    mean_x   = filter2(win, y, 'same');
    mean_y   = filter2(win, yn, 'same');
    cov_xy = filter2(win, y.*yn, 'same') - mean_x.*mean_y;
    ss_x  = filter2(win, y.^2, 'same') - mean_x.^2;
    ss_y  = filter2(win, yn.^2, 'same') - mean_y.^2;
    ss_x(ss_x<0)=0;
    ss_y(ss_y<0)=0;
    
    g = cov_xy./(ss_x+tol);  % estimate gain factor
    vv = (ss_y - g.*cov_xy); % estimation error
    g (ss_x < tol) = 0;
    vv (ss_x < tol) = ss_y (ss_x < tol);
    ss_x(ss_x<tol)=0;
    g (ss_y < tol) = 0;
    vv (ss_y < tol) = 0;
    
    aux = y;
    [Nsy,Nsx] = size(aux);
    prnt = parent & (nband < Nsc-1);   % check parent availability
    BL = zeros(size(aux,1),size(aux,2),1 + prnt);
    BL(:,:,1) = aux;
    if prnt,
        auxp = pyrBand(pyro, pind, nband+1);
        auxp = real(imenlarge2(auxp));
        BL(:,:,2) = auxp(1:Nsy,1:Nsx);
    end
    y=BL;
    [nv,nh,nb] = size(y);
    block = [blSzX blSzY];
      
    nblv = nv-block(1)+1;	% discard outer coefficients 
    nblh = nh-block(2)+1;   % for centrral coefficients (avoid boundary effects)
    nexp = nblv*nblh;		% number of coefficients
    N = prod(block) + prnt; % size of the neighborhood
    Ly = (block(1)-1)/2;
    Lx = (block(2)-1)/2;
    if (Ly~=floor(Ly))|(Lx~=floor(Lx)),
        error('Spatial dimensions of neighborhood must be odd!');
    end

    % Next: rearrange 'nexp' neighborhoods
    Y = zeros(nexp,N);
    n = 0;
    for ny=-Ly:Ly,	% spatial neighbors
        for nx=-Lx:Lx,
            n = n + 1;
            foo = shift(y(:,:,1),[ny nx]);
            foo = foo(Ly+1:Ly+nblv,Lx+1:Lx+nblh);
            Y(:,n) = foo(:);
        end
    end
    if prnt,	% parent
        n = n + 1;
        foo = y(:,:,2);
        foo = foo(Ly+1:Ly+nblv,Lx+1:Lx+nblh);
        Y(:,n) = foo(:);
    end

    C_u = innerProd(Y)/nexp;    % positive-definete covariance matrix
    [Q,L] = eig(C_u);           % eigenvalues with orthogonal matrix Q
    L = diag(diag(L).*(diag(L)>0))*sum(diag(L))/(sum(diag(L).*(diag(L)>0))+(sum(diag(L).*(diag(L)>0))==0)); % correct negative eigenvalues, maintaining variance
    C_u = Q*L*Q';
    ss = (Y*inv(C_u)).*Y/N;
    ss = sum(ss,2);
    ss = reshape(ss,nblv,nblh);
    L = diag(L);    
    g = g(Ly+1:Ly+nblv,Lx+1:Lx+nblh);
    vv = vv(Ly+1:Ly+nblv,Lx+1:Lx+nblh);
    
% compute info-weight using I(R;E|S)+I(R;F|S)-I(E;F|S), discarded
%     temp1=0; temp2=0; temp3=0;
%     for j=1:length(L)
%         temp1=temp1+(((log2(1+g.*g.*ss.*L(j)./(vv+sigma_nsq))))); % distorted image information for the i'th subband
%         temp2=temp2+(((log2(1+ss.*L(j)./(sigma_nsq))))); % reference image information
%  %       temp3=temp3+log2((ss.*L(j)+sigma_nsq).*(g.*g.*ss.*L(j)+vv+sigma_nsq)./(ss.*L(j).*(vv+sigma_nsq)+g.*g.*ss.*L(j)*sigma_nsq+sigma_nsq.*(sigma_nsq+vv))); % mutual information of E F
%        temp3=temp3+log2(1 + (ss.*g.*L(j)).^2./(ss.*(sigma_nsq+vv+g.*g.*sigma_nsq).*L(j) + sigma_nsq.*(vv+sigma_nsq))); % mutual information of E F
%     end
%     infow = temp1 + temp2 - temp3;

    infow = zeros(size(g));
    for j=1:length(L)
        infow = infow + log2(1 + ((vv+(1+g.*g).*sigma_nsq).*ss.*L(j)+sigma_nsq.*vv)./(sigma_nsq.*sigma_nsq)); % info-weight = I(R;E|S)+I(D;F|S)-I(E;F|S)
    end
    infow(infow < tol) = 0;
    iw_map{nband} = infow;
end

%%   
function imu = imenlarge2(im)

[M N] = size(im);
t1 = imresize(im, [4*M-3 4*N-3], 'bilinear');
t2 = zeros(4*M-1, 4*N-1);
t2(2:end-1, 2:end-1) = t1;
t2(1,:) = 2*t2(2,:) - t2(3,:);
t2(end,:) = 2*t2(end-1,:) - t2(end-2,:);
t2(:,1) = 2*t2(:,2) - t2(:,3);
t2(:,end) = 2*t2(:,end-1) - t2(:,end-2);
imu = t2(1:2:end, 1:2:end);