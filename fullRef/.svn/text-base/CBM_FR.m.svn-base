 
function FI=CBM_FR(origImg,distImg)
[areae areat areaf]=imseprate(origImg,distImg);

FSEe=FSEs_index(origImg,distImg,areae);

FSEt=FSEs_index(origImg,distImg,areat);
FSEf=FSEs_index(origImg,distImg,areaf);

ut=0.3375;
ue=0.4617;
uf=0.2008;

FI=(ue*FSEe+ut*FSEt+uf*FSEf)/(ue+ut+uf);

%%
function FSE=FSEs_index(origImg,distImg,area)

area=area(6:size(area,1)-5,6:size(area,2)-5);
   [s.mean s.map] = SSIM_FR(origImg,distImg);
   
   ssim=s.map(area>0);
   ssim=-sort(-ssim');
   N=size( ssim,2);
   
a0=1;
c0=size(ssim,2);
b0=round(a0+(c0-a0)/2);
As=ssim(a0);
Bs=ssim(b0);
Cs=ssim(c0);
a=a0;
b=b0;
c=c0;

Ag=size(find(ssim>=As),2)/N; 
Bg=size(find(ssim>=Bs),2)/N;
Cg=size(find(ssim>=Cs),2)/N;
while (c-a)>10
if Bs>Bg
        a=b;
        As=Bs;
        Ag=Bg;
    elseif Bs<Bg
        c=b;
        Cs=Bs;
        Cg=Bg;
    else
       break    
    end
    b=round(a+(c-a)/2);
    Bs=ssim(b);
    Bg=size(find(ssim>=Bs),2)/N;
end
if Bs==Bg
 FSE=Bs;
else

for i=a:c
    g(i)=size(find(ssim>=ssim(i)),2)/N;
    m(i)=min(g(i),ssim(i));
end
  FSE=max(m(i));
end

 
%%
function [areaE, areaT, areaF]=imseprate(imo,im)
%STAT Interesting statistics
% imo and im is the original and degraded  gray image respectively
% E=edge,F=flat,T=texture  
sobelh=[1 0 -1;
        2 0 -2;
        1 0 -1;];
sobelv=[-1 -2 -1;
         0 0 0;
         1 2 1 ];

imoh=filter2(sobelh,imo);
imov=filter2(sobelv,imo);
imoG=sqrt(imoh.^2+imov.^2);

imh=filter2(sobelh,im);
imv=filter2(sobelv,im);
imG=sqrt(imh.^2+imv.^2);


th1=0.16*max(max(imoG));
th2=0.08*max(max(imoG));

 [i,j,s] = find(imoG>th1);
 [m,n] = size(imoG);
 S1=sparse(i,j,s,m,n);
 clear i,clear j, clear s,clear m ,clear n
 [i,j,s] = find(imG>th1);
 [m,n] = size(imG);
 S2=sparse(i,j,s,m,n);
 areaE=full(S1+S2); 
 
  clear S1,clear S2;
 clear i,clear j, clear s,clear m ,clear n
 
 [i,j,s] = find(imoG<th2);
 [m,n] = size(imoG);
 S1 = sparse(i,j,s,m,n);
 clear i,clear j, clear s,clear m ,clear n
 [i,j,s] = find(imG<=th1);
 [m,n] = size(imG);
 S2=sparse(i,j,s,m,n);
areaF=full(S1.*S2);
  

  clear S1
 clear i,clear j, clear s,clear m ,clear n
   
 [i,j,s] = find((imoG>=th2)&(imoG<=th1));
 [m,n] = size(imoG);
 S1 = sparse(i,j,s,m,n);
 areaT=full(S1.*S2);
 