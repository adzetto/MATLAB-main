% Loudspeaker in closed cabinet model
npts=10000;
maxf=500;
f=maxf/npts:maxf/npts:maxf;			% freq range
w=2.*pi.*f;

% electrical
Re=6;			% 6 ohms
Le=0.003;		% 3 mH
Bl=13.25;				% Bl=6.5 gives Qt near 1, Bl=12 give Qt near .37

% diaphragm
fd=55;			% 25.5 Hz mech res
wd=2*pi*fd;
Md=0.060;			% 60 gm
Cd=1/(Md*wd^2);
Rd=1.5;

% cabinet
rho=1.21;
c=345;
Va=1;				% m^3 sealed box
Ca=Va/(rho*c^2);

% acoustic radiation
a=0.15;			% 12" is about 30cm
S=pi*a^2;
k=w./c;
i=sqrt(-1);
Za=(rho*c/S).*((.5.*k.*a).^2 + i.*0.6.*k.*a);
Ra=real(Za);
Ma=imag(Za);
Cad=Ca.*Cd./(S^2*Cd+Ca);

Zm=(Rd+(S^2).*Za) + i.*w.*Md + 1./(i.*w.*Cad);
Z1=Zm.*(Re+i.*w.*Le)+Bl^2;
Pae=(Bl.*S.*Ra)./Z1;

figure(1)
plot(f,20.*log10(abs(Pae)),'k');
xlabel('Hz');
ylabel('dB Pa / V');
title('Loudspeaker Transmitted Pressure Sensitivity');

% Theill=Small Params

w0sq=1/(Md*Cad);	% include cabinet
w0=sqrt(w0sq);
Qes=w0*Md*Re/(Bl)^2;
Qms=w0*Md/Rd;

Qt=w0*Md*Re/((Bl)^2+Re*Rd);	% total system Q

Zres=Re+Bl^2/Rd;					% max impedance
Vas=S^2*Cd*rho*c^2;
[Qes Qms Qt Vas*35.32 Zres]

% input electrical impedance

Z0=Re + i.*w.*Le + (i.*w.*Cd.*(Bl)^2)./(i.*w.*Cd.*Rd + (1. - w.^2./w0sq));

figure(2);

subplot(2,1,1),plot(f,(180/pi).*angle(Z0),'k');
title('Loudspeaker Input Electrical Impedance');
ylabel('Degrees');
subplot(2,1,2),plot(f,abs(Z0),'k');
ylabel('Ohms');
xlabel('Hz');

% microphone response

Ze=Bl^2./(Re+i.*w.*Le);
ep=S.*(Bl./(Ze + Rd + i.*w.*Md + 1./(i.*w.*Cd)));

figure(3)
plot(f,20.*log10(abs(ep)),'k');
xlabel('Hz');
ylabel('dB V / Pa');
title('Dynamic Microphone Voltage Sensitivity');
%axis([0 500 -60 0]);

% geophone response
gv=i.*w.*Md.*Bl./(Ze + Rd + i.*w.*Md + 1./(i.*w.*Cd));

figure(4)
plot(f,20.*log10(abs(gv)),'k');
xlabel('Hz');
ylabel('dB Vsec / m');
title('Geophone Voltage Sensitivity');


% loudspeaker impedance analysis

Zms=Rd + i.*w.*Md + 1./(i.*w.*Cd);
Zspk=Re + i.*w.*Le + Bl^2./Zms;

figure(5);
plot(f,abs(Zspk));

% estimate parameters from impedance curve

Le_est=(imag(Zspk(npts)-Zspk(round(.75*npts))))/(w(npts)-w(round(.75*npts)));
Re_est=abs(Zspk(1));
[Re_est Re Le_est Le]

[Zres,ires]=max(Zspk);
Zres=abs(Zres);

Zrsq=sqrt(Zres*Re_est);			% measure Qt
aZspk=abs(Zspk);
for n=2:ires,
   if aZspk(n) > Zrsq
      f1=f(n-1)+((aZspk(n)-aZspk(n-1))/(f(n)-f(n-1)))*(Zrsq-aZspk(n-1));
      break;
   end;
end;
for n=ires:npts,
   if aZspk(n) < Zrsq
      f2=f(n-1)+((aZspk(n)-aZspk(n-1))/(f(n)-f(n-1)))*(Zrsq-aZspk(n-1));
      break
   end;
end;
Qtest=(f(ires)/(f2-f1))*sqrt(Re/Zres);

% assume Md & Cd known from added-mass res shift
QtRd=w(ires)*Md*Re_est/Zres;
Rd_est=QtRd/Qtest;
Bl_est=sqrt((w(ires)*Md*Re_est/Qtest)-Re_est*Rd_est);

[Qt Qtest Bl Bl_est Rd Rd_est]


figure(6)

subplot(2,1,1),plot(f,(180/pi).*angle(Pae),'k');
title('Loudspeaker Transmitted Pressure Sensitivity');
ylabel('deg');
subplot(2,1,2),plot(f,abs(Pae),'k');
xlabel('Hz');
ylabel('Pa / V');
