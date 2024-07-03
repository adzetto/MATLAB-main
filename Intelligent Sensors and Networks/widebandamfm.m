% wideband AM and FM generator and detector

fs=20000;

fc=5000;
npts=8192;
n=1:npts;
t=n./fs;
f=0:fs/npts:fs-fs/npts;

wc=2*pi*fc/fs;

alpha=.9;
%beta=40;
beta=1;
beta=10;

fam=100;
ffm=40;
amsig=zeros(size(n));
fmsig=zeros(size(n));

% pulse train generator

% Nr/(2Nw) gives number of harmonics in freq domain!
Nr=400;	% sample # pulse repeats
Nw=20;		% pulse width in samples
tau=.5;	% small delay
% if beta is odd -> even harmonics
% if beta is even -> odd harmonics

Ts=1/fs;
t=Ts*n;

ypulse=(Nw/Nr).*sin(pi.*(n-tau)./Nw)./(sin(pi.*(n-tau)./Nr));

for nn=1:npts,
   if ypulse(nn)>1
      ypulse(nn)=1;
   end;
   if ypulse(nn)<-1
      ypulse(nn)=-1;
   end;
end;


y=(1. + amsig).*cos(wc.*n+beta.*ypulse);

Yf=sqrt(2).*abs(fft(y))./npts;

figure(1);
subplot(2,1,1),plot(t(1:400),y(1:400),'k');
title('FM Pulse Train');
ylabel('Response');
xlabel('Sec');
%axis([0 .01 -1 1]);
subplot(2,1,2),plot(f(1:4096),Yf(1:4096),'k');
ylabel('Response');
xlabel('Hz');

% apply Hilbert transform to new array

Yhf=zeros(size(Yf));

Yfc=fft(y);
Yhf=zeros(size(Yfc));
% bandpass filter and hilbert xform
nfc=2049;
band=2047;
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
npts=1500;
subplot(2,1,1),plot(t(1:npts),yfm(1:npts),'k');
title('FM Recovered Signal');
ylabel('Response');
xlabel('Sec');
%axis([0 .075 -1 1]);

subplot(2,1,2),plot(t(1:npts),yam(1:npts),'k');
title('AM Recovered Signal');
ylabel('Response');
xlabel('Sec');
axis([t(1) t(npts) -.1 .1]);
