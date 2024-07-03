% combo AM and FM generator and detector

fs=20000;

fc=5000;
fam=100;
ffm=40;

n=1:8192;
t=n./fs;
f=0:fs/8192:fs-fs/8192;

wc=2*pi*fc/fs;
wam=2*pi*fam/fs;
wfm=2*pi*ffm/fs;

alpha=.9;
beta=50;

y=(1. + alpha.*cos(wam.*n)).*cos(wc.*n + beta.*sin(wfm.*n));

Yf=sqrt(2).*abs(fft(y))./8192;

figure(1);
subplot(2,1,1),plot(t(1:300),y(1:300),'k');
title('FM=40Hz beta=50  AM=100Hz alpha=.9');
ylabel('Response');
xlabel('Sec');

subplot(2,1,2),plot(f(1000:3000),Yf(1000:3000),'k');
ylabel('Response');
xlabel('Hz');

% apply Hilbert transform to new array

Yhf=zeros(size(Yf));
band=200;	% 200 bins

Yfc=fft(y);
Yhf=zeros(size(Yfc));
% bandpass filter and hilbert xform
nfc=2049;
band=2000;
for n=-band:band,
   Yhf(nfc+n)=Yfc(nfc+n) + conj(Yfc(8192-nfc-n+2));  % 1 for matlab, 1 for dc
end;
yh=ifft(Yhf);

j=sqrt(-1);
n=1:8192;
ymod=exp(-j.*2.*pi.*fc.*n./fs).*yh;

yam=abs(ymod)-1.;
yfm=unwrap(angle(ymod));

figure(2);
npts=500;
subplot(2,1,1),plot(t(1:npts),yfm(1:npts),'k');
title('FM Recovered Signal');
ylabel('Response');
xlabel('Sec');

subplot(2,1,2),plot(t(1:npts),yam(1:npts),'k');
title('AM Recovered Signal');
ylabel('Response');
xlabel('Sec');
