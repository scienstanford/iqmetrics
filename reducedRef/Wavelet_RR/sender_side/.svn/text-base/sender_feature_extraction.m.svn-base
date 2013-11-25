
function NC1 = sender_feature_extraction(im,alpha)
                                 
 C=WTransform(double(im),[1,1,1]);%

n=0;


i=4;%(1)
w=WCSF(3/8);
% for j=1:2:8
    n=n+1;
    oband=C{i}{3};%(1)


    band=oband*w;
    th(n)=std2(band); 

thr=alpha*mean(th);

NC1(1)=thr;%NC1(1)


n=1;
i=2;
w=WCSF(3/32);
for j=1:3
%     for k= 1:2:8
      n=n+1;
      oband=C{i}{j};
      band=oband*w;
      absband=abs(band);
      NC1(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;
%     end                                                                   
end

i = 3;%(1)
% i=2;     %(2)
% n=1;     %(2)
w=WCSF(3/16);
for j=1:3
%   for k=2:2:8
     n=n+1;
     oband=C{i}{j};
     band=oband*w;
     absband=abs(band);
     NC1(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;
%   end
end
i = 4;%(1)

w=WCSF(3/8);
for j=1:3
%   for k=2:2:8
     n=n+1;
     oband=C{i}{j};
     band=oband*w;
     absband=abs(band);
     NC1(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;
%   end
end
return