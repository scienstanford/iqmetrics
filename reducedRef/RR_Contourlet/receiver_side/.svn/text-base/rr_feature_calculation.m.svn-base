function NC2 = rr_feature_calculation(im,NC1)
C=pdfbdec(double(im),'9-7','9-7',[3,3,4]);%下采样Contourlet
% C=nsctdec(double(im),[3,3,4]);%非下采样Contourlet
n=1;
thr=NC1(1);
i=2;
w=WCSF(3/32);
for j=1:2:8
    n=n+1;
    oband=C{i}{j};
    band=oband*w;
    absband=abs(band);
%     NC2(n) =round(100*length(find(absband>=thr))/prod(size(absband)))/100;
%     NC2(n) =round(1000*length(find(absband>=thr))/prod(size(absband)))/1000;
    NC2(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;
end

i=3;
w=WCSF(3/16);
for j=2:2:8
    n=n+1;
    oband=C{i}{j};
    band=oband*w; 
    absband=abs(band);
%     NC2(n) =round(100*length(find(absband>=thr))/prod(size(absband)))/100;
%     NC2(n) =round(1000*length(find(absband>=thr))/prod(size(absband)))/1000;
    NC2(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;
end

i = 4;
w=WCSF(3/8);
for j=1:2:16
    n=n+1;
    oband=C{i}{j};
    band=oband*w;
    absband=abs(band);
%     NC2(n) =round(100*length(find(absband>=thr))/prod(size(absband)))/100;
%     NC2(n) =round(1000*length(find(absband>=thr))/prod(size(absband)))/1000;
   NC2(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;
end
return