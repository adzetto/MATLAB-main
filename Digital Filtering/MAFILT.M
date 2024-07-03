% MA average filter

nwts=10;
B=(1/nwts).*ones([1 nwts]);

Zb=roots(B);
resp=fft(B,1024);
w=[1:512].*pi./512;

figure(1);

subplot(1,2,1), h=plot(real(Zb), imag(Zb),'ko');

axis([-1 1 -1 1]);
xlabel('Real z');
ylabel('Imag z');
title('10th Order MA Filter Zeros');

subplot(1,2,2), h=plot(w./(2*pi),abs(resp(1:512)),'k');

xlabel('Normalized Frequency f/fs');
ylabel('Linear Magnitude');
title('10 th Order MA Filter Frequency Response');

