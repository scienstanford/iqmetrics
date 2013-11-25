
function NC2 = receiver_distortion_measure(im,NC1)



 C=WTransform(double(im),[1,1,1]);  
%ãÐÖµ
thr=NC1(1);

n=1;
i=2;
w=WCSF(3/32);
for j=1:3
%   for k=1:2:8
      n=n+1;
      oband=C{i}{j};
      band=oband*w;
      absband=abs(band);
      NC2(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;
%    end
end

i = 3;%(1)
% i=2;    %(2)
% n=1;    %(2)
w=WCSF(3/16);
for j=1:3
%     for k=2:2:8
        n=n+1;
        oband=C{i}{j};
        band=oband*w;
        absband=abs(band);
        NC2(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;
%     end
end
i = 4;%(1)
% i=2;     %(2)
% n=1;     %(2)
w=WCSF(3/8);
for j=1:3
%   for k=2:2:8
     n=n+1;
     oband=C{i}{j};
     band=oband*w;
     absband=abs(band);
     NC2(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;
%   end
end
return