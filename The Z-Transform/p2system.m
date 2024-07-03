% response of simple conjugate pole pair

j=sqrt(-1);
fs=57;
T=1./fs;

fp1=25;
sigp1=-10;

% 1 conjg poles pair
sp1=sigp1+j.*2.*pi.*fp1;
sp1c=sigp1-j.*2.*pi.*fp1;

zp1=exp(sp1.*T);
zp1c=exp(sp1c.*T);

zscale=(zp1-zp1c)/(sp1-sp1c); % eqn 2.3.8
% impulse responses

As = 1 ./ (sp1-sp1c);
Az = 1 ./ (zp1-zp1c); 

n=1:300;
ta=n.*T;
hs=As.*exp(sp1.*n.*T)-As.*exp(sp1c.*n.*T);
hz=zscale.*( Az*(zp1.^(n-1))-Az.*(zp1c.^(n-1)) );

figure(1);
plot(ta,real(hz),'k.');
hold
plot(ta,real(hs),'k');
hold;

% system frequency responses

ifs=0:fs/200:fs/2;
s=0+j.*2.*pi.*ifs;
z=exp(s.*T);

hw=( 1./(sp1-sp1c) ).*( 1./(s-sp1) - 1./(s-sp1c) );  % eqn 2.2.13
hwzc= T.*( (zp1-zp1c)./(sp1-sp1c) ).*( 1./((z-zp1).*(z-zp1c)) ); % eqn 2.3.9
% note: multiply hwzc by z and they match, due to delay in digital system

figure(2);
subplot(2,1,1),plot(ifs,57.296.*(angle(hw)),'k',ifs,(180./pi).*(angle(hwzc)),'k.');
ylabel('Degrees');
xlabel('Hz');
legend('Equation 2.2.13','Equation 2.3.9');
title('Frequency Responses of Digital and Analog Systems');

subplot(2,1,2),plot(ifs,20.*log10(abs(hw)),'k',ifs,20.*log10(abs(hwzc)),'k.');
ylabel('dB');
xlabel('Hz');
legend('Equation 2.2.13','Equation 2.3.9');
