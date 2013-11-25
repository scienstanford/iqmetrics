function JND=JND_dct(distImg)

%%  DCT-based JND (Just-Noticeable Difference) model 

%This program implements the DCT-based JND (Just-Noticeable Difference) model 
%Xiaohui Zhang, Weisi Lin and Ping Xue, “Improved Estimation for Just-noticeable Visual Distortion”, 
%Signal Processing, Vol. 85(4), pp.795-808, April 2005.

%Related work to the pixel domain can be also found in:

%Xiaohui Zhang, Weisi Lin and Ping Xue, “Just-Noticeable Difference Estimation With Pixels in Images”, 
%Journal of Visual Communication and Image Representation, Vol 19(1), pp 30-41, 2008.
%Our proposed dct-jnd calculation, modified on 05/Jan/2004
%This program is to estimate the corresponding JND threshold values for
%each DCT coefficients of a given image, and then to generate a MND-noise
%contaminated image for evaluation. 
% I=imread('hats_JPEG.bmp');

% I = imread('lena.bmp');
% clear all;
% Image Read
% I = imread('frame1.tiff');
Lum = double(distImg);
k = 1;M =256;
[row,col] = size(distImg);
Lum=Lum(1:row/8*8,1:col/8*8);
% JND adjustment factor
tfac=0.3;

% Assume maximum luminance and minimum luminance values
Lmax = 130;
Lmin = 0;

%DCT Transform
Tr = dctmtx(8);
C = blkproc(Lum,[8 8],'P1*x*P2',Tr,Tr');

%conversion of Luminance-DCT coefficient
L = blkproc(C,[8 8], 'P1 + (P2 -P1)/P3*(x(1,1)/8 )',Lmin,Lmax,M);

%No. of 8X8 blocks in the image
[row1,col1] = size(L);

%% Parameters for CSF calculation
wx = 0.0298;  %horizontal pixel width
wy = 0.0298;   %vertical pixel width
r  = 0.6; LT = 13.45; S0 = 94.7; aT = 0.649; aL = 0.500; f0 = 6.78; af = 0.182; Lf = 300; K0 = 3.125; aK = 0.0706; LK = 300;
L0 = 65; s  = 0.25;
% Assume grey level 128 corresponding to mid-range luminance value 65 cd/m2
LB = 65;
C00 = 1024;
%Ahumada-Peterson CSF equation
if LB <= LK
    K = K0*((LB/LK)^aK);            
else 
    K = K0;
end
if LB <= Lf
    fmin = f0*((LB/Lf)^af);
else
    fmin = f0;
end
if LB <= LT
    Tmin = LT/S0*((LB/LT)^aT);
else
    Tmin = LB/S0;
end
for i = 1:8
    for j = 1:8
        freq(i,j) = sqrt(((i-1)^2)/(wx^2)+((j-1)^2)/(wy^2))/16;
        if freq(i,j) == 0
            sinang(1,1) = 0;
        elseif i == j
            sinang(i,j) = 1.0;
        else
            sinang(i,j) = 2*freq(i,1)*freq(1,j)/(freq(i,j)^2);
        end
    end
end
for i = 1:8
    for j = 1:8
        f = freq;
        f(1,1) = freq(1,2);
        ang(i,j) = asin(sinang(i,j));
        f1 = (log10(f(i,j)) - log10(fmin))^2;
        g(i,j) = log10(s*Tmin/(r+(1-r)*((cos(ang(i,j)))^2))) + K*f1;
        T(i,j) = 10.^(g(i,j));
    end
end
for i = 1:8
    for j = 1:8     
        if i == 1
            ai=sqrt(1/8);
        else
            ai=sqrt(2/8);
        end
        if j == 1
            aj=sqrt(1/8);
        else
            aj=sqrt(2/8);
        end
        tij(i,j) = M*T(i,j)/(Lmax-Lmin)/ai/aj;
    end
end
%% Adjust JND for DC
tij(1,1) = min(tij(1,2),tij(2,1));
C1=zeros(row1,col1,8,8);
for i = 1:row
    for j = 1:col
        C1(ceil(i/8),ceil(j/8),i+8-ceil(i/8)*8,j+8-ceil(j/8)*8) = C(i,j);
    end
end

%% Luminance Adaptation: a quasi-parabola function
aT=3; kT=2; aM=0.649; kQ=0.8; aQ=2;
for n1 = 1:row1
    for n2 = 1:col1
        for i=1:8
            for j=1:8
                if C1(n1,n2,1,1)>C00
                    tDCT(n1,n2,i,j) = tij(i,j)*(kQ*(C1(n1,n2)/C00-1)^aQ + 1);
                    aLum(n1,n2)=kQ*(C1(n1,n2)/C00-1)^aQ + 1;
                else
                    tDCT(n1,n2,i,j) = tij(i,j)*(kT*(1-C1(n1,n2)/C00)^aT + 1);
                    aLum(n1,n2)=kT*(1-C1(n1,n2)/C00)^aT + 1;
                end
            end
        end
    end
end

%% Block Classification
for n1 = 1:row1
    for n2 = 1:col1
        %texture masking model
        u1 = 125;
        u2 = 900;
        a1 = 2.3*3;
        b1 = 1.6*3;
        a2 = 1;
        b2 = 1.6;
        y1 = 2.0;
        y  = 4*4;
        k1  = 290;
        edge_area=1;
        texture_area=2;
        plain_area=0;
        edg(n1,n2) = sum(abs(C1(n1,n2,4:7,1)))+sum(abs(C1(n1,n2,1,4:7)))+sum(abs(C1(n1,n2,3,2:3)))...
                +abs(C1(n1,n2,2,3))+abs(C1(n1,n2,4,4));
        lowf(n1,n2) = sum(abs(C1(n1,n2,2:3,1)))+sum(abs(C1(n1,n2,1,2:3)))+abs(C1(n1,n2,2,2));
        highf(n1,n2) = sum(sum(abs(C1(n1,n2,:,:))))-edg(n1,n2)-lowf(n1,n2)-C1(n1,n2,1,1);
        
        edgn(n1,n2) = edg(n1,n2)/12;
        lowfn(n1,n2)=lowf(n1,n2)/5;
        highfn(n1,n2)=highf(n1,n2)/46;
        
        if edg(n1,n2) + highf(n1,n2) < u1
            block(n1,n2) = 'p';
        elseif edg(n1,n2) + highf(n1,n2) > u2
            if ((lowfn(n1,n2)/edgn(n1,n2) >= a2) & ((lowfn(n1,n2)+edgn(n1,n2))/highfn(n1,n2)>=b2))... 
                | ((lowfn(n1,n2)/edgn(n1,n2) >= b2) & ((lowfn(n1,n2)+edgn(n1,n2))/highfn(n1,n2)>=a2))...
                | ((lowfn(n1,n2)+edgn(n1,n2))/highfn(n1,n2)>=y1)
            block(n1,n2) = 'e';
            else block(n1,n2) = 't';
            end
        elseif ((lowfn(n1,n2)/edgn(n1,n2) >= a1) & ((lowfn(n1,n2)+edgn(n1,n2))/highfn(n1,n2)>=b1))...
                | ((lowfn(n1,n2)/edgn(n1,n2) >= b1) & ((lowfn(n1,n2)+edgn(n1,n2))/highfn(n1,n2)>=a1))...
                | ((lowfn(n1,n2)+edgn(n1,n2))/highfn(n1,n2)>=y)
            block(n1,n2) = 'e';
        elseif edg(n1,n2) + highf(n1,n2) > k1
            block(n1,n2) = 't';
        else block(n1,n2) = 'p';
        end        
        max1 = 1800;
        min1 = 290;
        if block(n1,n2) == 't'          %texture block
            TexE = edg(n1,n2) + highf(n1,n2);
%             FmaxT = 12;
            FmaxT = 2.5;
%             FmaxT = 2;
            TexMask(n1,n2) = (FmaxT -1)*(TexE - min1)/(max1 - min1) + 1;    
        elseif block(n1,n2) == 'e'      %edge block
            EdgE(n1,n2) = edg(n1,n2) + lowf(n1,n2);
            if EdgE(n1,n2) <= 400
                TexMask(n1,n2) = 1.125;
            else 
                TexMask(n1,n2) = 1.25;
            end
        elseif block(n1,n2) == 'p'      %plain block
            TexMask(n1,n2) = 1;
        end
    end
end

%% Texture masking elevation
for n1 = 1:row1
    for n2 = 1:col1
        for i = 1:8
            for j = 1:8
                if i+j == 2
                    aCM(n1,n2,1,1) = 1; %Masking elevation factor for DC
                else
                    aCM(n1,n2,i,j) = max(1,(abs(C1(n1,n2,i,j)/tDCT(n1,n2,i,j)))^0.36) * TexMask(n1,n2);
                end
                if (block(n1,n2) == 'e')  
                    aCM(n1,n2,2:7,1) = TexMask(n1,n2);
                    aCM(n1,n2,1,2:7) = TexMask(n1,n2);
                    aCM(n1,n2,2,2:3) = TexMask(n1,n2);
                    aCM(n1,n2,3,2:3) = TexMask(n1,n2);
                    aCM(n1,n2,4,4)   = TexMask(n1,n2);
                end
                %Final JND for each DCT coefficient         
                tJND(n1,n2,i,j) = tDCT(n1,n2,i,j)*aCM(n1,n2,i,j);
            end
        end
    end
end

    

for i = 1:row
    for j = 1:col
        JND(i,j) = tJND(ceil(i/8),ceil(j/8),i+8-ceil(i/8)*8,j+8-ceil(j/8)*8)*tfac;
    end
end
