% function y = WBCT(x,fname,dfilt,nlevs)
function y = WBCT(x,nlevs)
%Wavlet分解,x为输入图像，dfilt为分解滤波器组名,nlves为分解级数

if isempty(nlevs)
    y = {x};
else
    %%//////////////////////////////////////////////////////
    %首先进行小波分解
    x=double(x);
   [ca1,ch1,cv1,cd1]=dwt2(x,'db1');%%利用Daubechies滤波器
% %/////////////////////////////////////////////////////////////////////
%     [h0, h1] = dfilters('9-7', 'd');%小波变换用的滤波器名，（利用biorthogonal Daubechies 9-7滤波器）
%     [ca1,ch1,cv1,cd1]=wfb2dec(x,h0,h1);
    %%////////////////////////////////////////////////////////
    dfilt = '9-7';
     if nlevs(end) ~= 0
         % 对cd1进行nlevs(end)级方向滤波，得到2^nlevs(end)个方向子带
         %ch1后接方向滤波器组dfbdec
         %cv1后接方向滤波器组dfbdec
          switch dfilt        % Decide the method based on the filter name
            case {'pkva6', 'pkva8', 'pkva12', 'pkva'}   
                % Use the ladder structure (whihc is much more efficient)
                % 
                xcd_dir = dfbdec_l(cd1, dfilt, nlevs(end));
                xch_dir = dfbdec_l(ch1, dfilt, nlevs(end));
                xcv_dir = dfbdec_l(cv1, dfilt, nlevs(end));
            otherwise       
                % General case
                 xcd_dir = dfbdec(cd1, dfilt, nlevs(end)); %进行 nlevs层dfb分解,
                                                          %得到长为2^nlevs的单位向量上的子带图像    
                 xch_dir = dfbdec(ch1, dfilt, nlevs(end));
                 xcv_dir = dfbdec(cv1, dfilt, nlevs(end)); 
           end
           x_dir={xch_dir,xcv_dir,xcd_dir};
      end
       x_dir = {ch1, cv1, cd1};
       %然后对ca1部分反复进行上述分解
       ylo = WBCT(ca1,nlevs(1:end-1));
%         ylo = WBCT(ca1,fname, dfilt, nlevs(1:end-1));
       %最后得到分解结果
       y={ylo{:},x_dir};
end
   
       
         