% noise models

npts=1024;
fs=1024;
f=fs/npts:.5*fs;	% 0-512Hz
t=1/fs:1/fs:npts/fs;

% thermal& shot noise is white

xth=randn([1 npts]);
xfth=fft(xth)./npts;

figure(1);
subplot(2,1,1), plot(t,real(xth),'k');
ylabel('x[t]');
xlabel('Sec');
title('Thermal or Shot Noise');
subplot(2,1,2), plot(f,20.*log10(abs(xfth(1:.5*npts))),'k');
ylabel('dB');
xlabel('Hz');
axis([0 fs/2 -60 0]);

% contact noise variance is 1/f

xfr=zeros([1 npts]);
xfi=zeros([1 npts]);
xfc=zeros([1 npts]);
xt=zeros([1 npts]);

xfr(1:npts/2)=(1./f).*randn(1,.5*npts);
xfi(1:npts/2)=(1./f).*randn(1,.5*npts);

for k=1:.5*npts-1,
   xfr(npts-k+1)=xfr(k+1);
   xfi(npts-k+1)=-xfi(k+1);
end;
xfr(1)=1e-4;	% -80 dB
xfi(1)=0;

j=sqrt(-1);
xfc=xfr + j.*xfi;

xt=ifft(npts.*xfc);

figure(2);
subplot(2,1,1), plot(t,real(xt),'k');
ylabel('x[t]');
xlabel('Sec');
title('Contact Noise');
subplot(2,1,2), plot(f,20.*log10(abs(xfc(1:.5*npts))),'k');
ylabel('dB');
xlabel('Hz');
axis([0 fs/2 -80 0]);

% popcorn noise

xoff=zeros([1 npts]);

xoff(50:110)=5;	% pops...
xoff(230:350)=5;
xoff(410:450)=5;
xoff(700:760)=5;
xoff(780:870)=5;
xoff(930:975)=5;
xoff(990:1000)=5;

xtpop=.2.*xth + xoff;
xfpop=fft(xtpop)./npts;

figure(3);
subplot(2,1,1), plot(t,real(xtpop),'k');
ylabel('x[t]');
xlabel('Sec');
title('Popcorn Noise');
subplot(2,1,2), plot(f,20.*log10(abs(xfpop(1:.5*npts))),'k');
ylabel('dB');
xlabel('Hz');
axis([0 fs/2 -60 20]);

xoff(790:800)=5;
