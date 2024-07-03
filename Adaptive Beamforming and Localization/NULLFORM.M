% adaptive null former for linear array

nels=8;
nsnaps=10;

T=.1;	%time overlap for FFT's
f=100;
c=343;
lambda=c/f;
d=1.25;
w=2*pi*f;
k=w/c;
kd=w*d/c;

nsrc=3;
Ampsrc=[1 1 1];
Angsrc=[20 80 140];	%degrees
nz=1e-1;	% e-7 smallest, e+1 largest

srcang=Angsrc*pi/180.;

y0=zeros([nsnaps 1]);
Y=zeros([nsnaps nels-1]);
A=zeros([nels-1 1]);

j=sqrt(-1);

for n=1:nsnaps,
	phases=(2*pi).*([rand rand rand]);  % scr's random at each snapshot
	y0(n)=Ampsrc*exp(j*w*T*n)*exp(phases.') + nz*(randn+j*randn);
	for m=1:nels-1,
		Y(n,m)=Ampsrc*exp(j*w*T*n)*exp(-j.*kd*m.*cos(srcang.')+phases.') + nz*(randn+j*randn);
	end;
end;

A=inv(Y'*Y)*Y'*y0;
Ap=[1; -A];

z=roots(Ap);
angz=angle(z);
absz=abs(z);
theta=(180/pi).*acos(angz./kd);
[absz theta]
Angsrc

npts=1024;
AF=fftshift(fft(Ap,npts));
piscale=-pi:2*pi/npts:pi*(npts-1)/npts;
thetascale=(180/pi).*acos(piscale./kd);

figure(1);
plot(thetascale,20.*log10(abs(AF)),'k');
ax=axis;
%axis([0 180 -50 10]);
ylabel('dB');
xlabel('\theta');
