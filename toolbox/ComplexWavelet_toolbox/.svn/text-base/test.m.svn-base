x = imread('lena.bmp');
x = double(x);

J = 3;
[Faf, Fsf] = FSfarras;
[af, sf] = dualfilt1;
w0 = dualtree2D(x, J,Faf,af);
w1 = cplxdual2D(x, J, Faf, af);