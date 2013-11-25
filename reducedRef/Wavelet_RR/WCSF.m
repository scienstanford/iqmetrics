
function [Sw]=WCSF(w)
 l=0.8;
 v=61;
 fs=(2*l*v*tan(0.5*pi/180))/0.0254;
 f=w*fs;
 Sw=0.04992*(1+5.9375*f).*exp(-(0.114*f).^1.1);