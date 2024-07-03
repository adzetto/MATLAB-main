% compression/expansion amplifier responses

vin=[-10:.1:10];

a1=1;
a2=0.000000001;
a3=0.0035;

vlin=a1.*vin;
vexpand=a1.*vin + a2.*vin.^2 + a3.*vin.^3;
vcompress=a1.*vin + a2.*vin.^2 - a3.*vin.^3;

figure(1);
plot(vin,vlin,'k',vin,vexpand,'k--',vin,vcompress,'k-.');
ylabel('output');
xlabel('input');
legend('linear','expansion','compression','Location','SouthEast');
grid on;

% lets look at a sinusoid

n=[1:1024];
fs=16384;
f0=2000;
arg=2*pi*f0/fs;
y=10.*sin(arg.*n);
ycomp=a1.*y + a2.*y.^2 - a3.*y.^3;
yexpand=a1.*y + a2.*y.^2 + a3.*y.^3;

F=n.*fs./1024;
han=hanning(1024)';

Yf=fft(y.*han,1024);
Yfcomp=fft(ycomp.*han,1024);
Yfexpand=fft(yexpand.*han,1024);

figure(2);
plot(F(1:512),20.*log10(abs(Yf(1:512))));
title('linear');
figure(3);
plot(F(1:512),20.*log10(abs(Yfcomp(1:512))));
title('compress');
figure(4);
plot(F(1:512),20.*log10(abs(Yfexpand(1:512))));
title('expand');
