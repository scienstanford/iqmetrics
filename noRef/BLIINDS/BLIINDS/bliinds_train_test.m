%Michele Saad- BLIINDS- December 02,2009
%Code used to train and test BLIINDS on the entire LIVE database

a1=1; a2=1; a3=1; a4=1; a5=1; a6=1; a7=1; a8=1;
mtx=[min_kurt_all_L1(:,9).^a1 min_kurt_all_L2(:,9).^a2 min_cnt_all_L1(:,99).^a3 min_cnt_all_L2(:,99).^a4 mx_entropy_all_L1'.^a5 mx_entropy_all_L2'.^a6 variance_entropy_all_L1'.^a7 variance_entropy_all_L2'.^a8 dms_all];

clear count*;
count_train=0;
count_test=0;
for x=1:982
    if index_train_test(x)==1
        count_train=count_train+1;
        mtx_train(count_train,1)=mtx(x,1);
        mtx_train(count_train,2)=mtx(x,2);
        mtx_train(count_train,3)=mtx(x,3);
        mtx_train(count_train,4)=mtx(x,4);
        mtx_train(count_train,5)=mtx(x,5);
        mtx_train(count_train,6)=mtx(x,6);
        mtx_train(count_train,7)=mtx(x,7);
        mtx_train(count_train,8)=mtx(x,8);
        dms_all_train(count_train)=dms_all(x);
    else
        count_test=count_test+1;
        mtx_test(count_test,1)=mtx(x,1);
        mtx_test(count_test,2)=mtx(x,2);
        mtx_test(count_test,3)=mtx(x,3);
        mtx_test(count_test,4)=mtx(x,4);
        mtx_test(count_test,5)=mtx(x,5);
        mtx_test(count_test,6)=mtx(x,6);
        mtx_test(count_test,7)=mtx(x,7);
        mtx_test(count_test,8)=mtx(x,8);
        dms_all_test(count_test)=dms_all(x);
    end
end

a8=abs(corr(dms_all_test',mtx_test(:,8),'type','Spearman'));
a7=abs(corr(dms_all_test',mtx_test(:,7),'type','Spearman'));
a6=abs(corr(dms_all_test',mtx_test(:,6),'type','Spearman'));
a5=abs(corr(dms_all_test',mtx_test(:,5),'type','Spearman'));
a4=abs(corr(dms_all_test',mtx_test(:,4),'type','Spearman'));
a3=abs(corr(dms_all_test',mtx_test(:,3),'type','Spearman'));
a2=abs(corr(dms_all_test',mtx_test(:,2),'type','Spearman'));
a1=abs(corr(dms_all_test',mtx_test(:,1),'type','Spearman'));

sa=a1+a2+a3+a4+a5+a6+a7+a8;
 a1=a1/sa; a2=a2/sa; a3=a3/sa; a4=a4/sa; a5=a5/sa; a6=a6/sa; a7=a7/sa; a8=a8/sa;
   
 mtx=[min_kurt_all_L1(:,9).^a1 min_kurt_all_L2(:,9).^a2 min_cnt_all_L1(:,99).^a3 min_cnt_all_L2(:,99).^a4 mx_entropy_all_L1'.^a5 mx_entropy_all_L2'.^a6 variance_entropy_all_L1'.^a7 variance_entropy_all_L2'.^a8 dms_all];
 clear count*;
count_train=0;
count_test=0;
for x=1:982
    if index_train_test(x)==1
        count_train=count_train+1;
        mtx_train(count_train,1)=mtx(x,1);
        mtx_train(count_train,2)=mtx(x,2);
        mtx_train(count_train,3)=mtx(x,3);
        mtx_train(count_train,4)=mtx(x,4);
        mtx_train(count_train,5)=mtx(x,5);
        mtx_train(count_train,6)=mtx(x,6);
        mtx_train(count_train,7)=mtx(x,7);
        mtx_train(count_train,8)=mtx(x,8);
        dms_all_train(count_train)=dms_all(x);
    else
        count_test=count_test+1;
        mtx_test(count_test,1)=mtx(x,1);
        mtx_test(count_test,2)=mtx(x,2);
        mtx_test(count_test,3)=mtx(x,3);
        mtx_test(count_test,4)=mtx(x,4);
        mtx_test(count_test,5)=mtx(x,5);
        mtx_test(count_test,6)=mtx(x,6);
        mtx_test(count_test,7)=mtx(x,7);
        mtx_test(count_test,8)=mtx(x,8);
        dms_all_test(count_test)=dms_all(x);
    end
end

mu=mean([mtx_train dms_all_train']);
sig=cov([mtx_train dms_all_train']);

inc=0;
for j=1:10000
    dm_test(j)=inc;
    inc=inc+0.05;
end
 
for i=1:982-sum(index_train_test)      
        kt=[mtx_test(i,1).*ones(10000,1) mtx_test(i,2).*ones(10000,1) mtx_test(i,3).*ones(10000,1) mtx_test(i,4).*ones(10000,1) mtx_test(i,5).*ones(10000,1) mtx_test(i,6).*ones(10000,1) mtx_test(i,7).*ones(10000,1) mtx_test(i,8).*ones(10000,1)]; 
        pL(1:10000,i)=mvnpdf([kt dm_test'],mu,sig);
end
[tA,tB]=max(pL);
predictionN=dm_test(tB);

y=dms_all_test';
x=predictionN';

beta0 = [1/std(y), 1, mean(x), 1/std(y), mean(y)]; %ms
[beta_all,r,J,SIGMA,mse] = nlinfit(x,y,@fitting,beta0);
y2_all = fitting(beta_all,x);
all_corr = corr(y2_all,y)

scatter(dms_all_test,predictionN)
grid
title('Predicted DMOS VS Subjective DMOS','Fontsize',12)
xlabel('Subjective DMOS','Fontsize',12)
ylabel('Predicted DMOS','Fontsize',12)
