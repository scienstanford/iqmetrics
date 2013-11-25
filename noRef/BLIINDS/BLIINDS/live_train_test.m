%Michele Saad %December 02,2009
%Dividing LIVE according to CONTENT into TRAINING and TEST sets
clc;

load('/Users/michelesaad/Documents/Research/IQA DATABASES/LIVE-DATABASE/databaserelease2/refnames_all.mat');

train_names{1}='rapids.bmp';
train_names{2}='buildings.bmp';
train_names{3}='stream.bmp';
train_names{4}='churchandcapitol.bmp';
train_names{5}='lighthouse.bmp';
train_names{6}='coinsinfountain.bmp';
train_names{7}='lighthouse2.bmp';
train_names{8}='building2.bmp';
train_names{9}='paintedhouse.bmp';
train_names{10}='sailing1.bmp';
train_names{11}='sailing4.bmp';
train_names{12}='studentsculpture.bmp';
train_names{13}='flowersonih35.bmp';
train_names{14}='monarch.bmp';
train_names{15}='ocean.bmp';



test_names{1}='dancers.bmp';
test_names{2}='bikes.bmp';
test_names{3}='sailing3.bmp';
test_names{4}='carnivaldolls.bmp';
test_names{5}='woman.bmp';
test_names{6}='caps.bmp';
test_names{7}='statue.bmp';
test_names{8}='house.bmp';
test_names{9}='sailing2.bmp';
test_names{10}='cemetry.bmp';
test_names{11}='plane.bmp';
test_names{12}='parrots.bmp';
test_names{13}='womanhat.bmp';
test_names{14}='manfishing.bmp';
 
for i=1:982
    if strcmp(refnames_all{i},train_names{1}) || strcmp(refnames_all{i},train_names{2})|| strcmp(refnames_all{i},train_names{3}) || strcmp(refnames_all{i},train_names{4}) || strcmp(refnames_all{i},train_names{5}) || strcmp(refnames_all{i},train_names{6}) || strcmp(refnames_all{i},train_names{7}) || strcmp(refnames_all{i},train_names{8}) || strcmp(refnames_all{i},train_names{9}) || strcmp(refnames_all{i},train_names{10}) || strcmp(refnames_all{i},train_names{11}) || strcmp(refnames_all{i},train_names{12}) || strcmp(refnames_all{i},train_names{13}) || strcmp(refnames_all{i},train_names{14}) || strcmp(refnames_all{i},train_names{15}) %|| strcmp(refnames_all{i},train_names{16})   
        index_tt(i)=1;
    else
        index_tt(i)=0;
    end
end
% 
%  index_train_test=index_tt(1:227)';
%  index_train_test=index_tt(228:227+233)';
%  index_train_test=index_tt(1+227+233:174+227+233)';
%  index_train_test=index_tt(1+227+233+174:174+174+233+227)';
%  index_train_test=index_tt(1+227+233+174+174:end)';
   index_train_test=index_tt;
