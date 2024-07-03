% Rician pdf

npts=150;
rmax=10;
dr=rmax/npts;
r=0:dr:rmax-dr;

% signal levels (assume noise is unity variance)
A10=10;
A6=6;
A2=2;
A1=1;
A0=0;
oz=ones(size(r));

x0=0;
x1=A1.*r;
x2=A2.*r;
x6=A6.*r;
x10=A10.*r;

pm0=r.*exp(-(r.^2 + oz.*A0^2)./2).*besseli(0,x0);
p1=r.*exp(-(r.^2 + oz.*A1^2)./2).*besseli(0,x1);
p2=r.*exp(-(r.^2 + oz.*A2^2)./2).*besseli(0,x2);
p6=r.*exp(-(r.^2 + oz.*A6^2)./2).*besseli(0,x6);
p10=r.*exp(-(r.^2 + oz.*A10^2)./2).*besseli(0,x10);

figure(1);
nskip=10;   % skip and replot line for clear markings
plot(r(1:nskip:npts),pm0(1:nskip:npts),'ko',...
r(1:nskip:npts),p1(1:nskip:npts),'k+',...
r(1:nskip:npts),p2(1:nskip:npts),'k*',r,p6,'k');
legend('Noise Only','SNR 1','SNR 2','SNR 6');
ylabel('Probability Density');
xlabel('Received Waveform Magnitude');
title('Rician Densities for Various SNRs');
hold on;
plot(r,pm0,'k',r,p1,'k',r,p2,'k',r,p6,'k');
hold off

grid on;