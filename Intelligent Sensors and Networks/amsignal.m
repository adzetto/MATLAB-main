% AM modulation model

fs=20000;

fc=5000;
fsh=100;

n=1:8192;
t=n./fs;
f=0:fs/8192:fs-fs/8192;

wc=2*pi*fc/fs;
wsh=2*pi*fsh/fs;

y=(1. + 1.*cos(wsh.*n)).*cos(wc.*n);
Yf=sqrt(2).*abs(fft(y))./8192;

figure(1);
subplot(2,1,1),plot(t(1:200),y(1:200),'k');
title('AM Response');
ylabel('Response');
xlabel('Sec');

subplot(2,1,2),plot(f(1900:2200),Yf(1900:2200),'k');
ylabel('Response');
xlabel('Hz');

% lets look at multiple harmonic groups

fs=50000;
fc=5000;
fsh=350;

wc=2*pi*fc/fs;
wsh=2*pi*fsh/fs;
M=5;	% M sidebands
Q=4;	% Q mesh harmonics
fm=0:fs/8192:fs-fs/8192;

ym=zeros(size(y));
for q=1:Q,
   for m=1:M
      ym = ym + sqrt(1/q).*(1. + sqrt(1/m*q).*cos(m.*wsh.*n)).*cos(q.*wc.*n);
   end;
end;
Ymf=sqrt(2).*abs(fft(ym))./8192;

figure(2);
subplot(2,1,1),plot(t(1:200),ym(1:200),'k');
title('AM Response');
ylabel('Response');
xlabel('Sec');

subplot(2,1,2),plot(fm(1:4096),Ymf(1:4096),'k');
ylabel('Response');
xlabel('Hz');

