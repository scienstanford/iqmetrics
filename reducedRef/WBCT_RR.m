function quality = WBCT_RR(origImg,distImg)
path(path, 'toolbox/contourlet_toolbox');
alpha = 9;
NC1= sender_feature_extraction(origImg,alpha);

NC2 = receiver_distortion_measure(distImg,NC1);

D0=0.1;
k =sum(abs(NC1(2:end) - NC2(2:end)));

quality=1/(1+log(1+k/D0)/log(1.5));


%%
% function y = WBCT(x,fname,dfilt,nlevs)
function y = WBCT(x,nlevs)
%Wavlet·
if isempty(nlevs)
    y = {x};
else
 
    x=double(x);
   [ca1,ch1,cv1,cd1]=dwt2(x,'db1');%%
% %/////////////////////////////////////////////////////////////////////
%     [h0, h1] = dfilters('9-7', 'd');¡
%     [ca1,ch1,cv1,cd1]=wfb2dec(x,h0,h1);
    %%////////////////////////////////////////////////////////
    dfilt = '9-7';
     if nlevs(end) ~= 0

          switch dfilt        % Decide the method based on the filter name
            case {'pkva6', 'pkva8', 'pkva12', 'pkva'}   
                % Use the ladder structure (whihc is much more efficient)
                % 
                xcd_dir = dfbdec_l(cd1, dfilt, nlevs(end));
                xch_dir = dfbdec_l(ch1, dfilt, nlevs(end));
                xcv_dir = dfbdec_l(cv1, dfilt, nlevs(end));
            otherwise       
                % General case
                 xcd_dir = dfbdec(cd1, dfilt, nlevs(end)); 
                                                             
                 xch_dir = dfbdec(ch1, dfilt, nlevs(end));
                 xcv_dir = dfbdec(cv1, dfilt, nlevs(end)); 
           end
           x_dir={xch_dir,xcv_dir,xcd_dir};
      end
       x_dir = {ch1, cv1, cd1};
   
       ylo = WBCT(ca1,nlevs(1:end-1));

     
       y={ylo{:},x_dir};
end
   
       
         

%% CSF
function [Sw]=WCSF(w)
 l=0.8;
 v=61;
 fs=(2*l*v*tan(0.5*pi/180))/0.0254;
 f=w*fs;
 Sw=0.04992*(1+5.9375*f).*exp(-(0.114*f).^1.1);
 
%% sender side
function NC1 = sender_feature_extraction(origImg,alpha)

C=WBCT(origImg,[1,1,1]);
n=0;


i=4;%(1)
w=WCSF(3/8);
% for j=1:2:8
    n=n+1;
    oband=C{i}{3};%(1)


    band=oband*w;
    th(n)=std2(band);%Çó³ö±ê×¼²î 
% end
%/////////////////////////////////
%È·¶¨ãĞÖµ
% thr=1.5*mean(th);%%È·¶¨ãĞÖµ
% thr=mean(th);
% thr=2*mean(th);
% thr=3*mean(th);
% thr=4*mean(th);
thr=alpha*mean(th);
%//////////////////////////////
NC1(1)=thr;%NC1(1)ÎªãĞÖµ

%/////////////////////////////////////////////////////////////////
%(1)Ê×ÏÈ¼ÆËã´ÎÏ¸½Ú²ã£¨µÚ¶ş²ã£©ÖĞ¸ö·½ÏòÖĞÏµÊı´óÓÚãĞÖµµÄ±ÈÀıµÄ¡£
n=1;
i=2;
w=WCSF(3/32);%×Ó´ø¼ÓÈ¨ÏµÊı
for j=1:3
   for k= 1:2:8
      n=n+1;
      oband=C{i}{j};
      band=oband*w;
      absband=abs(band);
      NC1(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;%¸ø³öµÚ3²ã¸÷×Ó´ø´óÓÚãĞÖµËùÕ¼±ÈÀı
    end                                                                   % 
end
%/////////////////////////////////////////////////////////////////////////
%////////////////////////////////////////////////////////////////////////
%(1),(2)¼ÆËã×îÏ¸½Ú²ã£¨×îÍâ²ã£©ÖĞ¸ö·½ÏòÖĞÏµÊı´óÓÚãĞÖµµÄ±ÈÀıµÄ¡£
i = 3;%(1)
% i=2;     %(2)
% n=1;     %(2)
w=WCSF(3/16);
for j=1:3
  for k=2:2:8
     n=n+1;
     oband=C{i}{j};
     band=oband*w;
     absband=abs(band);
     NC1(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;%¸ø³öµÚ1²ãËùÕ¼±ÈÀı£¨×îÍâ²ã£©
  end
end
i = 4;%(1)

w=WCSF(3/8);
for j=1:3
   for k=2:2:8
     n=n+1;
     oband=C{i}{j};
     band=oband*w;
     absband=abs(band);
     NC1(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;%¸ø³öµÚ1²ãËùÕ¼±ÈÀı£¨×îÍâ²ã£©
  end
end


%%È·¶¨½µÖÊÍ¼ÏñÀï¸÷¸ö·½Ïò×Ó´øÖĞ·ùÖµ´óÓÚÒ»¶¨ãĞÖµÏµÊıÔÚ¸Ã×Ó´øÄÚËùÕ¼µÄ±ÈÀı
%NC1 ÓÉ²Î¿¼Í¼Ïñ»ñµÃ
%(1)WBCTÎª2²ã3¼¶·Ö½â
%(2)WBCTÎª1²ã3¼¶·Ö½â
function NC2 = receiver_distortion_measure(im,NC1)


%//////////////////////////////////////////////////////////////////
% C=WBCT(double(im),'db1','9-7',[3,3]);%(1)Ê×ÏÈ½øĞĞWBCT·Ö½â,Óë²Î¿¼¶Ë´¦ÀíÍêÈ«ÏàÍ¬
                               %ÎÒÃÇÈ¡Á½²ã·Ö½â 
 C=WBCT(double(im),[1,1,1]);  %(1)Ê×ÏÈ½øĞĞ3²ãĞ¡²¨·Ö½â,Óë²Î¿¼¶Ë´¦ÀíÍêÈ«ÏàÍ¬
                                  %ÎÒÃÇÈ¡Á½²ã·Ö½â                             
%  C=WBCT(double(im),'pkva',[3,3]);%(1)Ê×ÏÈ½øĞĞWBCT·Ö½â %ÎÒÃÇÈ¡Á½²ã·Ö½â                              
%/////////////////////////////////////////////////////////////////
% C=WBCT(double(im),'9-7',3);%(2)Ê×ÏÈ½øĞĞWBCT·Ö½â,Óë²Î¿¼¶Ë´¦ÀíÍêÈ«ÏàÍ¬
                               %ÎÒÃÇÈ¡1²ã3¼¶·Ö½â                                 

%////////////////////////////////////////////////////////////////
%ãĞÖµ
thr=NC1(1);
%////////////////////////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////
%(1)Óë²Î¿¼¶ËÏàÍ¬£¬Ê×ÏÈ¼ÆËã´ÎÏ¸½Ú²ã£¨µÚ¶ş²ã£©ÖĞ¸ö·½ÏòÖĞÏµÊı´óÓÚãĞÖµµÄ±ÈÀıµÄ¡£
n=1;
i=2;
w=WCSF(3/32);
for j=1:3
  for k=1:2:8
      n=n+1;
      oband=C{i}{j};
      band=oband*w;
      absband=abs(band);
      NC2(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;
   end
end
%/////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////
%(1),(2)Óë²Î¿¼¶ËÏàÍ¬£¬¼ÆËã×îÏ¸½Ú²ã£¨µÚÈı²ã£©ÖĞ¸ö·½ÏòÖĞÏµÊı´óÓÚãĞÖµµÄ±ÈÀıµÄ¡£
i = 3;%(1)
% i=2;    %(2)
% n=1;    %(2)
w=WCSF(3/16);
for j=1:3
    for k=2:2:8
        n=n+1;
        oband=C{i}{j};
        band=oband*w;
        absband=abs(band);
        NC2(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;
    end
end
i = 4;%(1)
% i=2;     %(2)
% n=1;     %(2)
w=WCSF(3/8);
for j=1:3
  for k=2:2:8
     n=n+1;
     oband=C{i}{j};
     band=oband*w;
     absband=abs(band);
     NC2(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;%¸ø³öµÚ1²ãËùÕ¼±ÈÀı£¨×îÍâ²ã£©
  end
end
return