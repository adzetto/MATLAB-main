% FM signal generator

fs=20000;
fc=5000;
fm=100;

beta=.5;

n=1:8196;
t=n./fs;
f=0:fs/8192:fs-fs/8192;

wc=2*pi*fc/fs;
wm=2*pi*fm/fs;

y=cos(wc.*n + beta.*sin(wm.*n));

Yf=sqrt(2).*abs(fft(y))./8192;

figure(1);
subplot(2,1,1),plot(t(1:200),y(1:200),'k');
title('FM Response beta=.5');
ylabel('Response');
xlabel('Sec');

subplot(2,1,2),plot(f(1900:2200),Yf(1900:2200),'k');
ylabel('Response');
xlabel('Hz');

beta=40;
y=cos(wc.*n + beta.*sin(wm.*n));

Yf=sqrt(2).*abs(fft(y))./8192;

figure(2);
subplot(2,1,1),plot(t(1:200),y(1:200),'k');
title('FM Response beta=40');
ylabel('Response');
xlabel('Sec');

subplot(2,1,2),plot(f(1:4096),Yf(1:4096),'k');
ylabel('Response');
xlabel('Hz');
