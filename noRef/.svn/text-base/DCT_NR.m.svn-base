function proportion = DCT_NR(distImg)


 W=[   0.4247    0.7326    0.9150    0.9789    0.9621    0.8964    0.8047    0.7029
       0.7326    0.8518    0.9505    0.9808    0.9511    0.8813    0.7896    0.6893
       0.9150    0.9505    0.9789    0.9713    0.9252    0.8507    0.7598    0.6628
       0.9789    0.9808    0.9713    0.9387    0.8813    0.8047    0.7168    0.6250
       0.9621    0.9511    0.9252    0.8813    0.8200    0.7452    0.6628    0.5782
       0.8964    0.8813    0.8507    0.8047    0.7452    0.6759    0.6011    0.5251
       0.8047    0.7896    0.7598    0.7168    0.6628    0.6011    0.5352    0.4686
       0.7029    0.6893    0.6628    0.6250    0.5782    0.5251    0.4686    0.4114];
% 
% W1 = 1./W;
% W1 = W1/sum(sum(W1))*100;
% W1 = W/sum(sum(W))*100;
W(1,1) = 0;  % don't consider dc coefficient 

[row,cols] = size(distImg);

blkwidth = 8;
row_blk = floor(row/blkwidth);
cols_blk = floor(cols/blkwidth);

% check whether the image size is integer divided by block width or not
if(blkwidth*row_blk~=row || blkwidth*cols_blk~=cols)
    img = distImg(1:blkwidth*row_blk,1:blkwidth*cols_blk);
%     img_ref =img_ref(1:blkwidth*row_blk,1:blkwidth*cols_blk);
    [row,cols] = size(img);
end

% S = blkproc(img_ref,[blkwidth blkwidth],@std2);

% Y = blkproc(img,[blkwidth blkwidth],@dct2);
Y = blkproc(distImg-128,[blkwidth blkwidth],@dct2);

cont = zeros(blkwidth,blkwidth);
% cont = zeros(row_blk,cols_blk);
gamma = 81;
labda = log(1+(255-gamma)^0.5)/log(1+(gamma)^0.5);

for ii = 1:row_blk
    for jj = 1:cols_blk
        tmp = Y((ii-1)*blkwidth+1:ii*blkwidth,(jj-1)*blkwidth+1:jj*blkwidth);
%         tt = img((ii-1)*blkwidth+1:ii*blkwidth,(jj-1)*blkwidth+1:jj*blkwidth);
%         %         if(S(ii,jj)<th1)
%         %             cont1 = cont1 + sum(sum(abs(tmp)<1e-5));
%         %         elseif(S(ii,jj)>th1 && S(ii,jj)<th2)
%         %             cont2 = cont2 + sum(sum(abs(tmp)<1e-5));
%         %         else
%         %             cont3 = cont3 + sum(sum(abs(tmp)<1e-5));
%         %         end
%         cont(ii,jj) = sum(sum(abs(tmp)<1e-05));
%         dc = tmp(1,1);
%         if(abs(dc)<1e-05)
%             cont(ii,jj) = cont(ii,jj)-1;
%         end
%         
%         lum = max([dc/8+128,0]);
%         lum = min([255,lum]);
%         
%         deta = std2(tt);
%         
%         if(lum<gamma)
%             wi = labda *log(1+lum^0.5/(1+deta));
%         else
%             wi= log(1+(255-lum)^0.5/(1+deta));
%         end
%         cont(ii,jj) = cont(ii,jj)*wi;
        for m = 1:blkwidth
            for n = 1:blkwidth
                cont(m,n) = cont(m,n) + (abs(tmp(m,n))<1e-05);
            end
        end
    end
end

% wi = wi/sum(sum(wi)); 
% cont = cont.*wi;

% proportion = (w1*cont1+w2*cont2+w3*cont3)/(row*cols);
proportion = sum(sum((cont/(row*cols)).*W));
