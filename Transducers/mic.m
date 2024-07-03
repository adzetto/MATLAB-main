% condenser mics

f=logspace(0,5,200);
w = 2.*pi.*f;
rho = 1.21;
c = 345;
k = w./c;

% electrical

Sd	= 4e-6;	% 4mm^2
a = sqrt(Sd/pi);
d	= 1e-4;	% .1 mm
Vb	= 36;	% 4 9V batteries
phi = Vb/d;	% V/m or N/C
e0	=	8.85e-12;
C0	= Sd*e0/d;
R0 = 1e30;

% diaphragm

Md	= 1e-3;	% 1 gram
f0 = 20000;
Cd	= 1/(4*pi^2*Md*f0^2);
Rd	= 1e2;

% Cavity
Vcav = .01*.01*.005;	% (1cm)^3
Cav = Vcav/(Sd^2*rho*c^2);

% Vent
Lv=0.001;	% 1mm long
rv=0.00005;	% .05 mm radius
eta=181e-6;	% air viscosity in Poise
Zv=8*eta*Lv*Sd^2/(pi*rv^4) + i.*w.*(4*rho*Lv*Sd^2/(3*pi*rv^2));

Zmr = (rho*c/Sd).*i.*.6.*k.*a;	% acs mass loading

Zm = 1./(i.*w.*Cd) + i.*w.*Md + Rd + Zmr.*Sd^2;
Ze=((C0*phi)^2).*R0./(1 + i.*w.*R0.*C0);

%ep = (1000*C0*phi)./((i.*w.*C0 + 1/R0).*( 1./(i.*w.*Cav) + Zm ) +(C0*phi)^2);

ep=1000.*Ze.*Zv./((C0*phi./(i.*w.*Cav)).*(Zv+Zm+Ze) + (C0*phi)*Zv.*(Zm+Ze));

rad2deg=180/pi;

% vary vent diameter
figure(1);
subplot(2,1,1),semilogx(f,rad2deg.*angle(ep),'k');
title('Condenser Mic Responses vs Vent Diameter');
ylabel('deg');
hold on;
subplot(2,1,2),semilogx(f,abs(ep),'k');
hold on;

rv=0.0001;	% .1 mm radius
Zv=8*eta*Lv*Sd^2/(pi*rv^4) + i.*w.*(4*rho*Lv*Sd^2/(3*pi*rv^2));
ep=1000.*Ze.*Zv./((C0*phi./(i.*w.*Cav)).*(Zv+Zm+Ze) + (C0*phi)*Zv.*(Zm+Ze));
subplot(2,1,1),semilogx(f,rad2deg.*angle(ep),'k:');
subplot(2,1,2),semilogx(f,abs(ep),'k:');

rv=0.0002;	% .2 mm radius
Zv=8*eta*Lv*Sd^2/(pi*rv^4) + i.*w.*(4*rho*Lv*Sd^2/(3*pi*rv^2));
ep=1000.*Ze.*Zv./((C0*phi./(i.*w.*Cav)).*(Zv+Zm+Ze) + (C0*phi)*Zv.*(Zm+Ze));
subplot(2,1,1),semilogx(f,rad2deg.*angle(ep),'k--');
axis([1 1e5 -200 200]);
legend('.1mm vent','.2mm vent','.4mm vent');
hold off;
subplot(2,1,2),semilogx(f,abs(ep),'k--');
hold off;

ylabel('mV/Pa');
xlabel('Hz');

% vary R0

figure(2);
rv=0.0000001;	% no radius
Zv=8*eta*Lv*Sd^2/(pi*rv^4) + i.*w.*(4*rho*Lv*Sd^2/(3*pi*rv^2));

R0=1e12;	% 1000 gig
Ze=((C0*phi)^2).*R0./(1 + i.*w.*R0.*C0);
ep=1000.*Ze.*Zv./((C0*phi./(i.*w.*Cav)).*(Zv+Zm+Ze) + (C0*phi)*Zv.*(Zm+Ze));

subplot(2,1,1),semilogx(f,rad2deg.*angle(ep),'k');
title('Condenser Mic Responses vs Shunt Resistance R_0');
ylabel('deg');
hold on;
subplot(2,1,2),semilogx(f,abs(ep),'k');
hold on;

R0=1e11;	% 100G
Ze=((C0*phi)^2).*R0./(1 + i.*w.*R0.*C0);
ep=1000.*Ze.*Zv./((C0*phi./(i.*w.*Cav)).*(Zv+Zm+Ze) + (C0*phi)*Zv.*(Zm+Ze));
subplot(2,1,1),semilogx(f,rad2deg.*angle(ep),'k:');
subplot(2,1,2),semilogx(f,abs(ep),'k:');

R0=1e10;	% 10g
Ze=((C0*phi)^2).*R0./(1 + i.*w.*R0.*C0);
ep=1000.*Ze.*Zv./((C0*phi./(i.*w.*Cav)).*(Zv+Zm+Ze) + (C0*phi)*Zv.*(Zm+Ze));
subplot(2,1,1),semilogx(f,rad2deg.*angle(ep),'k--');
axis([1 1e5 -200 200]);
legend('1 Tera-ohm','100 Giga-ohms','10 Giga-ohms');
hold off;
subplot(2,1,2),semilogx(f,abs(ep),'k--');
hold off;

ylabel('mV/Pa');
xlabel('Hz');