% basic Fourier resolution

fs=100;
npts=256;
n=1:npts;

%f0=fs/4;
f0=28.125;


y=sin((2*pi*f0/fs).*(n-1));  
% must start @ 0 for real-imag to match theory when f0=fs/4
dzp=16; % play with this to see different resolutions

T=dzp./fs;  % total seconds of integration time
df=1/T;
titbuf=sprintf('T=%.3f df=%.3f',T,df);

for k=dzp+1:npts,    % zero pad
    y(k)=0;
end;

NT=dzp;
Yf=(1/NT).*fftshift(fft(y,npts));
Yfdzp=(1/NT).*fftshift(fft(y(1:dzp),dzp));

F=(n-npts/2).*fs./npts;
Fdzp=((1:(npts/dzp):npts)-npts/2).*fs./npts;

Fn=F./fs;
Fndzp=Fdzp./fs;

figure(1);
plot(Fn,abs(Yf),'k',Fndzp,abs(Yfdzp),'kO');
xlabel('f/f_s');
ylabel('|Y(f)/T|');
title(titbuf);
axis([-.5 .5 0 .6]);

figure(2);
plot(Fn,real(Yf),'k',Fndzp,real(Yfdzp),'kO');
xlabel('Hz');
ylabel('Re{Y(f)/T}');
title(titbuf);
axis([-.5 .5 -.6 .6]);

figure(3);
plot(Fn,imag(Yf),'k',Fndzp,imag(Yfdzp),'kO');
%plot(F./f0,imag(Yf),'k');   % show curve only for Figure 1
xlabel('f/f_s');
%xlabel('f/f0');
ylabel('imag{Y(f)/T}');
title(titbuf);
axis([-.5 .5 -.6 .6]);

