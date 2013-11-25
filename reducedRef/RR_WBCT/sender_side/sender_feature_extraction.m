%确定原图像里各个方向子带中幅值大于一定阈值系数在该子带内所占的比例，并给出阈值NC1: thr
% function NC1 = feature_calculation(im)
function NC1 = sender_feature_extraction(im,alpha)
%alpha为阈值加权值
%(1)WBCT 2层3级分解
%(2)WBCT 1层3级分解

%/////////////////////////////////////////////////////////////
% C=WBCT(double(im),'db1','9-7',[3,3]);%(1)首先进行WBCT分解 ，小波方向滤波都用9-7
                                       %我们取两层分解
C=WBCT(im,[1,1,1]);%%(1)首先进行3层小波分解 ，小波方向滤波都用9-7
                                       %我们取两层分解
%//////////////////////////////////////////////////////////////
% C=WBCT(double(im),'9-7',3);%(2)首先进行WBCT分解 %我们取1层分解
%/////////////////////////////////////////////////////////////
n=0;

%//////////////////////////////
%确定阈值
i=4;%(1)
w=WCSF(3/8);
% for j=1:2:8
    n=n+1;
    oband=C{i}{3};%(1)


    band=oband*w;
    th(n)=std2(band);%求出标准差 
% end
%/////////////////////////////////
%确定阈值
% thr=1.5*mean(th);%%确定阈值
% thr=mean(th);
% thr=2*mean(th);
% thr=3*mean(th);
% thr=4*mean(th);
thr=alpha*mean(th);
%//////////////////////////////
NC1(1)=thr;%NC1(1)为阈值

%/////////////////////////////////////////////////////////////////
%(1)首先计算次细节层（第二层）中个方向中系数大于阈值的比例的。
n=1;
i=2;
w=WCSF(3/32);%子带加权系数
for j=1:3
   for k= 1:2:8
      n=n+1;
      oband=C{i}{j};
      band=oband*w;
      absband=abs(band);
      NC1(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;%给出第3层各子带大于阈值所占比例
    end                                                                   % 
end
%/////////////////////////////////////////////////////////////////////////
%////////////////////////////////////////////////////////////////////////
%(1),(2)计算最细节层（最外层）中个方向中系数大于阈值的比例的。
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
     NC1(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;%给出第1层所占比例（最外层）
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
     NC1(n) =round(512*length(find(absband>=thr))/prod(size(absband)))/512;%给出第1层所占比例（最外层）
  end
end