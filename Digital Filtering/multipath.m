% multipath model for Figure 3.2

fs=10000;
T=1/fs;
npts=1024;
n=1:npts;
f=(fs/npts).*n;

x=zeros([1 npts]);
x2=x;

x(1)=1;
x(15)=.95;

x2(1)=1;
x2(15)=.5;

x3(1)=1;
x3(15)=.1;

Xf=fft(x,1024);
Xf2=fft(x2,1024);
Xf3=fft(x3,1024);

figure(1);
plot(f(1:npts/2),20.*log10(abs(Xf(1:npts/2))),'k',...
    f(1:npts/2),20.*log10(abs(Xf2(1:npts/2))),'k.',...
    f(1:npts/2),20.*log10(abs(Xf3(1:npts/2))),'k:');
ylabel('dB');
xlabel('Hz');
axis([0 fs/2 -25 10]);
legend('95% reflection','50% reflection','10% reflection');