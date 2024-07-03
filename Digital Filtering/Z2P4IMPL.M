%
% analog-digital system simulation Mar 19, 1994 DCS
% 
j=sqrt(-1);

fs=1000;		% analog frequency limit
T=1./fs;		% time step size
fd=600;		% digital sample rate
Td=1./fd;		% digital sample interval
ntimepts=200;		% samples in impulse response plots
nfreqpts=200;		% samples in frequency response plots

% 1 conjg zeros pair, 2 conjg pole pairs
fz1=130;
sigz1=+5;
fp1=240;
sigp1=-20;
fp2=160;
sigp2=-10;

sz1=sigz1+j.*2.*pi.*fz1;	% s-plane poles and zeros
sz1c=sigz1-j.*2.*pi.*fz1;
sp1=sigp1+j.*2.*pi.*fp1;
sp1c=sigp1-j.*2.*pi.*fp1;
sp2=sigp2+j.*2.*pi.*fp2;
sp2c=sigp2-j.*2.*pi.*fp2;

%
% Partial fraction expansion coeffs analog system
%
As=((sp1-sz1)*(sp1-sz1c))/((sp1-sp1c)*(sp1-sp2)*(sp1-sp2c));
Bs=((sp1c-sz1)*(sp1c-sz1c))/((sp1c-sp1)*(sp1c-sp2)*(sp1c-sp2c));
Cs=((sp2-sz1)*(sp2-sz1c))/((sp2-sp1)*(sp2-sp1c)*(sp2-sp2c));
Ds=((sp2c-sz1)*(sp2c-sz1c))/((sp2c-sp1)*(sp2c-sp1c)*(sp2c-sp2));

% analog system impulse response
n=0:ntimepts-1;
h=As*exp(sp1.*n.*Td)+Bs*exp(sp1c.*n.*Td)+Cs*exp(sp2.*n.*Td)+Ds*exp(sp2c.*n.*Td);

% analog system frequency response
ifs=0:fd/(2.*nfreqpts):fd/2;
s=0+j.*2.*pi.*ifs;
hw=((s-sz1).*(s-sz1c))./((s-sp1).*(s-sp1c).*(s-sp2).*(s-sp2c));

%
% map s-plane system to digital domain z-plane
%

% z-plane poles and zeros
zz1=exp(sz1.*Td);	
zz1c=exp(sz1c.*Td);
zp1=exp(sp1.*Td);
zp1c=exp(sp1c.*Td);
zp2=exp(sp2.*Td);
zp2c=exp(sp2c.*Td);

%
% Partial fraction expansion coeffs digital z-plane system
%
Az=((zp1-zz1)*(zp1-zz1c))/((zp1-zp1c)*(zp1-zp2)*(zp1-zp2c));
Bz=((zp1c-zz1)*(zp1c-zz1c))/((zp1c-zp1)*(zp1c-zp2)*(zp1c-zp2c));
Cz=((zp2-zz1)*(zp2-zz1c))/((zp2-zp1)*(zp2-zp1c)*(zp2-zp2c));
Dz=((zp2c-zz1)*(zp2c-zz1c))/((zp2c-zp1)*(zp2c-zp1c)*(zp2c-zp2));

%
% Partial fraction expansion coeffs for z-plane all-pole system
%
Ac=1./((sp1-sp1c)*(sp1-sp2)*(sp1-sp2c));
Bc=1./((sp1c-sp1)*(sp1c-sp2)*(sp1c-sp2c));
Cc=1./((sp2-sp1)*(sp2-sp1c)*(sp2-sp2c));
Dc=1./((sp2c-sp1)*(sp2c-sp1c)*(sp2c-sp2));

% digital system impulse response
k=1:ntimepts;
tdd=(k+1).*Td;
ta=n.*Td;

% mapped amplitude scaled digital impulse response
hz=Td.*(Az*zp1.^k+Bz*zp1c.^k+Cz*zp2.^k+Dz*zp2c.^k);
%hz=exp(sigz1.*Td).*(1./Td^2).*(Ac*zp1.^k+Bc*zp1c.^k+Cc*zp2.^k+Dc*zp2c.^k);
% mapped impulse response with modal scaling
hzc=As*zp1.^k+Bs*zp1c.^k+Cs*zp2.^k+Ds*zp2c.^k;
%
% plot impulse responses
%
figure(1);
% unscaled
plot(tdd,real(hz),'r.');	
hold on
% analog system
plot(ta,real(h),'y');		
% modal scaled response
plot(tdd,real(hzc),'g.');	
hold off;
%
% write to disk
%
tresp = [ta; real(h); real(hz); real(hzc)];
fd = fopen('z2p4t6.dat','wt');
fprintf(fd,'%e\t%e\t%e\t%e\n',tresp);
fclose(fd);

%
% Digital Frequency Responses
%
z=exp(s.*Td);			

% mapped amplitude scaled digital response
hwz=Td.*((z-zz1).*(z-zz1c)).*(1./Td^2).*((Ac./(z-zp1))+(Bc./(z-zp1c))+(Cc./(z-zp2))+(Dc./(z-zp2c)));

% mapped modal scaled digital response
hwzc=Td.*((As./(z-zp1))+(Bs./(z-zp1c))+(Cs./(z-zp2))+(Ds./(z-zp2c)));
%
% plot frequency responses
%
figure(2);
% amplitude scaled
plot(ifs,20.*log10(abs(hwz)),'r');	
hold on
% analog response
plot(ifs,20.*log10(abs(hw)),'y');	
% modal scaled
plot(ifs,20.*log10(abs(hwzc)),'g');	
hold off
%
% write to disk
%
fresp = [ifs; (20.*log10(abs(hw)));(20.*log10(abs(hwz)));(20.*log10(abs(hwzc)))];
fd = fopen('z2p4h6.dat','wt');
fprintf(fd,'%e\t%e\t%e\t%e\n',fresp);
fclose(fd);




