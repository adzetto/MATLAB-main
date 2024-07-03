% PZT model

d11 = 1000e-12;					% 1000 pC/N typical PZT
phi = 1/d11;

% stud compliance and mass
din=0.125;						% stud diameter in inches
S=pi*(din*0.0254*.5)^2;		% x-section area in m^2
lin=0.25;
l=lin*0.0254;
rho=7750;						% stainless density & speed
c=7000;

Cs=l/(rho*S*c^2);
Ms=rho*l*S;

% PZT compliance and mass
din=0.25;						% PZT diameter in inches
S=pi*(din*0.0254*.5)^2;		% x-section area in m^2
lin=0.25;
l=lin*0.0254;
rho=1750;						% PZT density & speed
c=1e7;

Cc=l/(rho*S*c^2);
Mc=rho*l*S;

Mp = 0.0025;						% 2.5 gr proof mass
R=1000;							% relatively small mech damping

f=logspace(2,10,200);
w=2.*pi.*f;
wsq=w.^2;

% electrical impedance
C0=1600e-12;							% 1600 pF electrical cap

num = (1./(i.*w.*Cc) + (C0*phi^2)./(i.*w)).*(i.*w.*(.5*Mc+.3*Ms+Mp) + 1./(i.*w.*Cs)+R);
den = 1./(i.*w.*Cc) + (C0*phi^2)./(i.*w) + i.*w.*(.5*Mc+.3*Ms+Mp) + 1./(i.*w.*Cs)+R;
Z1= num./den;
Zm=i.*w.*Mc.*0.5 + Z1;
Mt=Mc+Ms+Mp;
u=Mt./Zm;
f1=u.*Z1;
um=f./(1./(i.*w.*Cc) + (C0*phi^2)./(i.*w));
ech=phi.*u./(i.*w);
ea=9810.*ech;	% mV/g

figure(1)
semilogx(f,abs(ea),'k');
ylabel('mV / g');
xlabel('Hz');
title('Typical Accelerometer Voltage Sensitivity');
axis([100 1000000 0 90]);
% hydrophone

rhow = 64*.45359/.028317;		% 1 lb = .45359 Kg, 1 ft^3 is .028317 m^3
cw=1500;
a = .5*din.*.0254;
S=pi*a^2;
k = w./cw;
Zam=(rhow*cw/S^3).*((.5.*k.*a).^2 + i.*.6.*k.*a);
Ram=real(Zam);

Zm=1./(i.*w.*Cc)+i.*w.*.5.*Mc + R;
Zet=1./(i.*w.*C0) + (phi^2./w.^2)./(Zm + Zam);
Zl=i.*w;	% 1 Kg load
Zetl=1./(i.*w.*C0) + (phi^2./w.^2)./(Zm + Zl);

pe = .001.*Ram./(S.*Zet.*(C0*phi + (i.*w./phi).*(Zm + Zam)));
f_i = .001.*Zl./(C0*phi + (i.*w./phi).*(Zm + Zl));

figure(2);
semilogx(f,abs(pe),'k');
ylabel('Pa/mV');
xlabel('Hz');
title('Hydrophone Transmited Pressure Response')

figure(3);
semilogx(f,20.*log10(abs(f_i)),'k');
ylabel('dB re 1 Nt/mA');
xlabel('Hz');
title('Force Actuator Response - 1 Kg Load');

% force gauge or receiver

ep=-1000.*phi./(i.*w.*S.*((C0*phi^2)./(i.*w) + Zm));

figure(4);
semilogx(f,20.*log10(abs(ep)),'k');
ylabel('dB re mV/Pa');
xlabel('Hz');
title('PZT Hydrophone Receiving Sensitivity');