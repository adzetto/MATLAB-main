% dipole velocity movie

rc=415;

%f=input('freq?');
f=171.5;

%f=1000;
%d=input('d?');
d=1;

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

% source 1 at < -.5d, 0 > source 2 at < +.5d, 0 >

R1=sqrt((X+.5*d*ones(size(X))).^2+Y.^2+1.e-10);		% ranges to field points from 1
R2=sqrt((X-.5*d*ones(size(X))).^2+Y.^2+1.e-10);		% ranges to field points from 2

figure(1);

nframes=16
fs=nframes*f;
T=1/fs;

% initial velocity calculation

TH1=atan2(Y,X+.5*d*ones(size(X)));	% angles to field points from source 1
TH2=atan2(Y,X-.5*d*ones(size(X)));	% angles to field points from source 2

A=j.*k.*rc.*Rscl;
Q=Q.*Rscl;

n=0;

U1=A.*exp(-j.*k.*R1).*(ones(size(R1))-j./(k.*R1))./(rc.*R1);
U1X=cos(TH1).*U1;
U1Y=sin(TH1).*U1;

U2=A.*exp(-j.*k.*R2).*(ones(size(R2))-j./(k.*R2))./(rc.*R2);
U2X=cos(TH2).*U2;
U2Y=sin(TH2).*U2;

UX=U1X - U2X;
UY=U1Y - U2Y;

% pressure calculation
P1 =j.*k.*rc.*Q.*exp(-j.*k.*R1)./(4.*pi.*R1);
P2 =-j.*k.*rc.*Q.*exp(-j.*k.*R2)./(4.*pi.*R2);
P = P1 + P2;

% instantaneous intensity
IX=UX.*P;
IY=UY.*P;

figure(1);
quiver(X,Y,real(IX),real(IY),'k');
lim=axis;

M=moviein(nframes);

for n=1:nframes,

% velocity frames

 U1=A.*exp(j*2*pi*f*n*T).*exp(-j.*k.*R1).*(ones(size(R1))-j./(k.*R1))./(rc.*R1);
 U1X=cos(TH1).*U1;
 U1Y=sin(TH1).*U1;

 U2=A.*exp(j*2*pi*f*n*T).*exp(-j.*k.*R2).*(ones(size(R2))-j./(k.*R2))./(rc.*R2);
 U2X=cos(TH2).*U2;
 U2Y=sin(TH2).*U2;

 UX=U1X - U2X;
 UY=U1Y - U2Y;

 % pressure calculation
P1 =j.*k.*rc.*Q.*exp(-j.*k.*R1)./(4.*pi.*R1);
P2 =-j.*k.*rc.*Q.*exp(-j.*k.*R2)./(4.*pi.*R2);
P = P1 + P2;

% instantaneous intensity
IX=UX.*P;
IY=UY.*P;
 
 quiver(X,Y,real(IX),real(IY),'k');
 axis(lim);

 M(:,n)=getframe;

end
movie(M,3);
	
