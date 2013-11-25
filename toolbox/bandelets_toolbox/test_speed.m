% test for haar transform speed
%
%   Copyright (c) 2005 Gabriel Peyré

n = 1024;
p = 1000;
x = rand(n,p);

tic;
y1 = perform_haar_transform_slow(x,1);
toc;

tic;
y2 = zeros(n,p);
for i=1:p
  y2(:,i) = perform_haar_transform(x(:,i),1);    
end
toc;