% nearfield intensity

kr=logspace(0,3,25);

krs=2*pi*.1;	% source radius is 1/10 of a wavelength
krm=2*pi;	% source radius is 1 wavelengths
krl=2*pi*10000;	% source radius is 10 wavelengths

j=sqrt(-1);

u0=0.01;	% 10 cm/sec
rc=415;		% rho-c is 415 Rayls

Aks=j*(krs^2)*u0*rc*exp(j*krs)/(1+j*krs);
Akm=j*(krm^2)*u0*rc*exp(j*krm)/(1+j*krm);
Akl=j*(krl^2)*u0*rc*exp(j*krl)/(1+j*krl);

Ints=abs(Aks^2).*(ones(size(kr)) + j./kr)./(2.*(kr.^2).*rc);
Intm=abs(Akm^2).*(ones(size(kr)) + j./kr)./(2.*(kr.^2).*rc);
Intl=abs(Akl^2).*(ones(size(kr)) + j./kr)./(2.*(kr.^2).*rc);


figure(1);
subplot(2,1,1), semilogx(kr,180.*unwrap(angle(Ints))./pi,'ko');
hold on;
subplot(2,1,1), semilogx(kr,180.*unwrap(angle(Intm))./pi,'k+');
subplot(2,1,1), semilogx(kr,180.*unwrap(angle(Intl))./pi,'k');
hold off;
ylabel('Degrees');
legend('\lambda/10','1 \lambda','10 \lambda');
title('Intensity Phase and Power vs Source Size and Distance');

subplot(2,1,2), loglog(kr,abs(Ints),'ko');
hold on;
subplot(2,1,2), loglog(kr,abs(Intm),'k+');
subplot(2,1,2), loglog(kr,abs(Intl),'k');
hold off;

% nearfield does not depend on source size!

xlabel('Distance in Wavelengths');
ylabel('Watts/sq meter');

% compute acoustic impedance

az=rc.*1./(ones(size(kr))+1./(j.*kr));

figure(2);
subplot(2,1,1), semilogx(kr,180.*unwrap(angle(az))./pi,'k');
ylabel('Degrees');
title('Radiation Impedance of Spherical Wave');
subplot(2,1,2), semilogx(kr,abs(az),'k');
xlabel('Distance in Wavelengths');
ylabel('Rayls');



R=1:1000;
lambda=1;
dphi=(180/pi).*2.*pi.*(R./lambda).*(1-cos(atan(lambda./(2.*R))));

figure(3);
semilogx(R,dphi,'k');
axis([10 1000 0 5]);
ylabel('Degrees');
xlabel(' Distance in R/\lambda from Point Source');
title('Phase Error Over 1-\lambda patch from Plane Wave Assumption');