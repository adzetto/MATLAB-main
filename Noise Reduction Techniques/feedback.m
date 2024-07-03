% 1 HG Feedback analysis

npts=1000;
f=logspace(0,5,npts);

w=(2*pi).*f;

fc=300;

wc=2*pi*fc;

Rf=100;	% unity gain filter resistor
Cf=1/(Rf*wc);

j=sqrt(-1);

Gf=1./(ones(size(f)) + j.*w.*Cf.*Rf);

% 2, 4, and 8-pole filters

Gf2=Gf.*Gf;
Gf4=Gf2.*Gf2;
Gf8=Gf4.*Gf4;

% note subtraction assumed at summing junction

figure(1);
subplot(2,1,1),semilogx(f,(180/pi).*angle(Gf),'k');
ylabel('deg');
title('300 Hz Low-Pass Filter Responses');
hold on;
subplot(2,1,1),semilogx(f,(180/pi).*angle(Gf2),'k:');
subplot(2,1,1),semilogx(f,(180/pi).*angle(Gf4),'k--');
subplot(2,1,1),semilogx(f,(180/pi).*angle(Gf8),'k-.');
hold off;
%legend('1-pole','2-pole','4-pole','8-pole');

subplot(2,1,2),semilogx(f,20.*log10(abs((Gf))),'k');
ylabel('dB');
xlabel('Hz');
hold on;
subplot(2,1,2),semilogx(f,20.*log10(abs((Gf2))),'k:');
subplot(2,1,2),semilogx(f,20.*log10(abs((Gf4))),'k--');
subplot(2,1,2),semilogx(f,20.*log10(abs((Gf8))),'k-.');
hold off;
legend('1-pole','2-pole','4-pole','8-pole');
axis([f(1) f(npts) -100 0]);


% y/n is 1/(1 + HS)

%tau=0.0;
%gain=1;
tau=input('delay:');
gain=input('gain:');

H=gain.*exp(-j.*w.*tau);

Yn1=1./(ones(size(f))+H.*Gf);
Yn2=1./(ones(size(f))+H.*Gf2);
Yn4=1./(ones(size(f))+H.*Gf4);
Yn8=1./(ones(size(f))+H.*Gf8);

figure(2);
subplot(2,1,1),semilogx(f,(180/pi).*angle(Yn1),'k');
ylabel('deg');
tbuff=sprintf('Closed-Loop Response G=%.f T=%.3f ms',gain, tau*1000);
title(tbuff);
hold on;
subplot(2,1,1),semilogx(f,(180/pi).*angle(Yn2),'k:');
subplot(2,1,1),semilogx(f,(180/pi).*angle(Yn4),'k--');
subplot(2,1,1),semilogx(f,(180/pi).*angle(Yn8),'k-.');
hold off;
v=axis;
axis([v(1) v(2) -181 181]);
%legend('1-pole','2-pole','4-pole','8-pole');

subplot(2,1,2),semilogx(f,20.*log10(abs((Yn1))),'k');
ylabel('dB');
xlabel('Hz');
hold on;
subplot(2,1,2),semilogx(f,20.*log10(abs((Yn2))),'k:');
subplot(2,1,2),semilogx(f,20.*log10(abs((Yn4))),'k--');
subplot(2,1,2),semilogx(f,20.*log10(abs((Yn8))),'k-.');
hold off;
legend('1-pole','2-pole','4-pole','8-pole');

% look at HG

HG1=H.*Gf;
HG2=H.*Gf2;
HG4=H.*Gf4;
HG8=H.*Gf8;

figure(3);
plot(real(HG1),imag(HG1),'k');
title('Open Loop H[s]G[s] Response');
hold on;
plot(real(HG2),imag(HG2),'k:');
plot(real(HG4),imag(HG4),'k--');
plot(real(HG8),imag(HG8),'k-.');
plot(-1,0,'k*');
hold off;
ylabel('imag [HG]');
xlabel('real [HG]');
axis([-gain gain -gain gain]);
legend('1-pole','2-pole','4-pole','8-pole','[-1,0]');
