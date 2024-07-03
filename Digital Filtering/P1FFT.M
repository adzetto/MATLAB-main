fs=1024;
T=1/fs;
npts=1024;
n=0:T:(npts-1)*T;
f=0:fs/npts:fs*(1-1/npts);

fp1=374;
sigp1=-10;

j=sqrt(-1);

sp1=sigp1+j.*2.*pi.*fp1;
sp1c=sigp1-j.*2.*pi.*fp1;

%
% digital domin
%

zp1=exp(sp1.*T);
zp1c=exp(sp1c.*T);

A=zeros(1024,1);
A(1)=1.0;
A(2)=-(zp1+zp1c);
A(3)=zp1.*zp1c;

AW=fft(A);
HAW=1./AW;
%
% digital filter
%
X=randn(1024,1);
Y=zeros(1024,1);
H=zeros(1024,1);

Y(1)=A(1).*X(1);
Y(2)=A(1).*X(2)-A(2).*Y(1);
for i=3:1024,
	Y(i)=A(1).*X(i) - A(2).*Y(i-1) - A(3).*Y(i-2);
end

YW=fft(Y);
XW=fft(X);
HW=XW./YW;

H=ifft(HW);
for i=6:1024,
	H(i)=0.0;
end

HW=fft(H);

HA=ifft(HAW);
HA32=HA;
HA32(33:1024)=zeros(size(HA32(33:1024)));
HAW32=fft(HA32);
HA8=HA;
HA8(9:1024)=zeros(size(HA8(9:1024)));
HAW8=fft(HA8);

figure(1);
plot(f,20.*log10(abs(HAW)),'k');
hold on;
plot(f,20.*log10(abs(HAW32)),'k.');
plot(f,20.*log10(abs(HAW8)),'k-.');
hold off;
axis([0 fs/2 -40 40]);
legend('1024 points','32 points', '8 points');
xlabel('Hz');
ylabel('dB');
