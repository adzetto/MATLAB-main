% dipole velocity movie

rc=415;

%f=input('freq?');
f=171.5;

%f=1000;
%d=input('d?');
d=1;
c=343;
%r0=.0127;	% 1" approx diameter sources
r0=.1;
d0=.0005;	% 1mm p-p displacement
u0=2*pi*f*d0;	% surface velocity at given f

k=2*pi*f/343;
kd=k*d;
j=sqrt(-1);

Q=4*pi*(r0^2)*u0*(1-sin(kd)/kd);	% dipole pressure coupling

max=3;		% grid range and spacing
step=d/4;

[X,Y]=meshgrid(-max-step/2:step:max+step/2,-max-step/2:step:max+step/2);
Rscl=sqrt(X.^2+Y.^2);  % This is to scale the arrow lengths
TH=atan2(Y,X);	% angles to field points from origin

A=j.*k.*rc;
Q=Q;

R1=sqrt((X+.5*d*ones(size(X))).^2+Y.^2+1.e-10);		% ranges to field points from 1
R2=sqrt((X-.5*d*ones(size(X))).^2+Y.^2+1.e-10);		% ranges to field points from 2
TH1=atan2(Y,X+.5*d*ones(size(X)));	% angles to field points from source 1
TH2=atan2(Y,X-.5*d*ones(size(X)));	% angles to field points from source 2

% pressure fields from each source
P1 =j.*k.*rc.*Q.*exp(-j.*k.*R1)./(4.*pi.*R1);
P2 =j.*k.*rc.*Q.*exp(-j.*k.*R2)./(4.*pi.*R2);

% velocity field of each source
U1=A.*exp(-j.*k.*R1).*(ones(size(R1))-j./(k.*R1))./(rc.*R1);
U1X=U1.*cos(TH1);
U1Y=U1.*sin(TH1);

I1=.5.*P1.*conj(U1);
I1log=log(abs(I1./1e-12)).*I1./abs(I1);    %preserve sign make size log scale
I1logX=I1log.*cos(TH1);                     % 1e-12 is 0 dB
I1logY=I1log.*sin(TH1);

U2=A.*exp(-j.*k.*R2).*(ones(size(R2))-j./(k.*R2))./(rc.*R2);
U2X=U2.*cos(TH2);
U2Y=U2.*sin(TH2);

I2=.5.*P2.*conj(U2);
I2X=I2.*cos(TH2);
I2Y=I2.*sin(TH2);
I2log=log(abs(I2./1e-12)).*I2./abs(I2);    %preserve sign make size log scale
I2logX=I2log.*cos(TH2);                     % 1e-12 is 0 dB
I2logY=I2log.*sin(TH2);

figure(1);
quiver(X,Y,real(I1logX),real(I1logY),'k');
figure(2);
quiver(X,Y,real(I2logX),real(I2logY),'k');

figure(3);
PT=P1-P2;
UTX=U1X+U2X;
UTY=U1Y+U2Y;
ITX=.5.*PT.*conj(UTX);
ITY=.5.*PT.*conj(UTY);

ITlogX=log(abs(ITX./1e-12)).*ITX./abs(ITX);    %preserve sign make size log scale
ITlogY=log(abs(ITY./1e-12)).*ITY./abs(ITY);
quiver(X,Y,real(ITlogX),real(ITlogY),'k');

% what a pain!
% now lets just plot the time averaged intensity correctly

I=(k^2*c*Q^2/(8*pi^2)).*(1./Rscl).*(ones(size(Rscl))-j./(k.*Rscl)).*sin(.5.*k.*d.*cos(TH)).^2;
Ilog=log(abs(real(I)./1e-12)).*real(I)./abs(I);
IlogX=Ilog.*cos(TH);
IlogY=Ilog.*sin(TH);
IX=I.*cos(TH);
IY=I.*sin(TH);
figure(4);
quiver(X,Y,real(IX),real(IY),'k');

