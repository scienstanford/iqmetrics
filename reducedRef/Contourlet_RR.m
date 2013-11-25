function quality = Contourlet_RR(origImg,distImg)

path(path, 'toolbox/contourlet_toolbox');
alpha = 9;
NC1= sender_feature_extraction(origImg,alpha);

NC2 = receiver_distortion_measure(distImg,NC1);

D0=0.1;
k =sum(abs(NC1(2:end) - NC2(2:end)));

quality=1/(1+log(1+k/D0)/log(1.5));


%% sender_feature_extraction
function NC1 = sender_feature_extraction(origImg,alpha)
                                 
C=pdfbdec(double(origImg),'9-7','9-7',[3,3,4]);%

n=0;

i=4;
w=WCSF(3/8);
for j=1:2:16
    n=n+1;
    oband=C{i}{j};
    band=oband*w;
 
    th(n)=std2(band);
end

thr=alpha*mean(th);

NC1(1)=thr;

n=1;
i=2;
w=WCSF(3/32);
for j=1:2:8
   
    n=n+1;
    oband=C{i}{j};
    band=oband*w;
    absband=abs(band);

    NC1(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;
end

i=3;
w=WCSF(3/16);
for j=2:2:8
    n=n+1;
    oband=C{i}{j};
    band=oband*w; 
    absband=abs(band);

    NC1(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;
end

i = 4;
w=WCSF(3/8);
for j=1:2:16
    n=n+1;
    oband=C{i}{j};
    band=oband*w;
    absband=abs(band);

    NC1(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;
end
return

%% CSF
function [Sw]=WCSF(w)
 l=0.8;
 v=61;
 fs=(2*l*v*tan(0.5*pi/180))/0.0254;
 f=w*fs;
 Sw=0.04992*(1+5.9375*f).*exp(-(0.114*f).^1.1);
 
%%
function NC2 = receiver_distortion_measure(distImg,NC1)
C=pdfbdec(double(distImg),'9-7','9-7',[3,3,4]);
n=1;
thr=NC1(1);
i=2;
w=WCSF(3/32);
for j=1:2:8
    n=n+1;
    oband=C{i}{j};
    band=oband*w;
    absband=abs(band);

    NC2(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;
end

i=3;
w=WCSF(3/16);
for j=2:2:8
    n=n+1;
    oband=C{i}{j};
    band=oband*w; 
    absband=abs(band);

    NC2(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;
end

i = 4;
w=WCSF(3/8);
for j=1:2:16
    n=n+1;
    oband=C{i}{j};
    band=oband*w;
    absband=abs(band);

   NC2(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;
end
return