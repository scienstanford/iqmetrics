function [q1,q2]=FFTJpeg(im)
%input:
%       im:gray of the JPEG image
%
%output:
%       q:quality of the JPEG image
%       p_v:vertical
%       p_h:horizontal
%
%
%-----------------code begin-----------------------------
im=double(im);
N=512;
K=4;
%vertical
g=abs(im(:,2:end)-im(:,1:end-1));

g=g';
s=reshape(g,1,size(g,1)*size(g,2));

n=size(s,2);
L=floor(n/N);

for ii=1:L
    ss=s((ii-1)*N+1:ii*N);
    ff_s=abs(fft(ss));

    p_v(ii,:)=2*(ff_s(1:N/2+1).^2);
    p_v(ii,1)=ff_s(1).^2;
    p_v(ii,N/2+1)=ff_s(N/2+1).^2;
end

p_v=sum(p_v)/L;
p_v=log10(1+p_v);
p_v=p_v/max(p_v);
p_v1=p_v;
pp=zeros(1,4);
iii=1;
for ii=1+K:N/2-K
    s_pv=sort(p_v1(ii-K:ii+K));
    if(ii==N/8+1 | ii==2*N/8+1 | ii==3*N/8+1 | ii==4*N/8+1) %论文中是N是从0开始
        pp(iii)=s_pv(K+1);
        iii=iii+1;
        continue;
    end
    p_v(ii)=s_pv(K+1);
end

M_bv=(p_v(N/8+1)-pp(1)+p_v(2*N/8+1)-pp(2)+p_v(3*N/8+1)-pp(3)+p_v(4*N/8+1)-pp(4))*8/7;
 
%horizontal
g=abs(im(2:end,:)-im(1:end-1,:));

s=reshape(g,1,size(g,1)*size(g,2));

n=size(s,2);
L=floor(n/N);

for ii=1:L
    ss=s((ii-1)*N+1:ii*N);
    ff_s=abs(fft(ss));

    p_h(ii,:)=2*(ff_s(1:N/2+1).^2);
    p_h(ii,1)=ff_s(1).^2;
    p_h(ii,N/2+1)=ff_s(N/2+1).^2;
end

p_h=sum(p_h)/L;
p_h=log10(1+p_h);
p_h=p_h/max(p_h);
p_h1=p_h;
pp=zeros(1,4);
iii=1;
for ii=1+K:N/2-K
    s_ph=sort(p_h1(ii-K:ii+K));
    if(ii==N/8+1 | ii==2*N/8+1 | ii==3*N/8+1 | ii==4*N/8+1) %论文中是N是从0开始
        pp(iii)=s_ph(K+1);
        iii=iii+1;
        continue;
    end
    p_h(ii)=s_ph(K+1);
end

M_bh=(p_h(N/8+1)-pp(1)+p_h(2*N/8+1)-pp(2)+p_h(3*N/8+1)-pp(3)+p_h(4*N/8+1)-pp(4))*8/7;

%all
 q1=(M_bv+M_bh)/2;
 q2=max(M_bv,M_bh);
 
%----------------the end---------------------