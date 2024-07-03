% Doppler wavelet

CPA=200;
v=100;

Ttotal=20;
f0=2;
c=344;
fs=1000;
nhalf=round(.5*Ttotal*fs);
n=1:2*nhalf+1;
T=(n-nhalf)./fs;
j=sqrt(-1);

r=sqrt(CPA^2+v^2.*T.^2);
f=f0.*( 1-(2/c).*(v^2).*T./r );
h=(1./r).*exp((j*2*pi).*f.*T);

figure(1);
titbuf=sprintf('CPA=%.f F0=%.f V=%.f',CPA,f0,v);
subplot(2,1,1),plot(T,real(h));
title(titbuf);
ylabel('h[t]');
xlabel('t sec');

CPA=20;
v=100;
r=sqrt(CPA^2+v^2.*T.^2);
f=f0.*( 1-(2/c).*(v^2).*T./r );
h=(1./r).*exp((j*2*pi).*f.*T);
titbuf=sprintf('CPA=%.f F0=%.f V=%.f',CPA,f0,v);
subplot(2,1,2),plot(T,real(h));
title(titbuf);
ylabel('h[t]');
xlabel('t sec');
