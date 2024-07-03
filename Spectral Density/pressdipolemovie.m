% dipole pressure movie

rc=415;

%f=input('freq?');
f=171.5;
f=343;
%d=input('d?');
d=1;

titbuf=sprintf('d=%.2f f=%.2f',d,f);

r0=.0127;	% 1" approx diameter sources
d0=.0005;	% 1mm p-p displacement
u0=2*pi*f*d0;	% surface velocity at given f

k=2*pi*f/343;
kd=k*d;
j=sqrt(-1);

Q=4*pi*(r0^2)*u0*(1-sin(kd)/kd);	% dipole pressure coupling

max=3;		% grid range and spacing
step=.155;

[X,Y]=meshgrid(-max:step:max,-max:step:max);

% source 1 at < -.5d, 0 > source 2 at < +.5d, 0 >

R1=sqrt((X+.5*d*ones(size(X))).^2+Y.^2+1.e-10);		% ranges to field points from 1
R2=sqrt((X-.5*d*ones(size(X))).^2+Y.^2+1.e-10);		% ranges to field points from 2

figure(1);
colormap(bone);

nframes=64;
fs=nframes*f;
T=1/fs;

% initial pressure calculation

P1 =j.*k.*rc.*Q.*exp(-j.*k.*R1)./(4.*pi.*R1);
P2 =-j.*k.*rc.*Q.*exp(-j.*k.*R2)./(4.*pi.*R2);
P = P1 + P2;
surfl(X,Y,imag(P));
shading interp;
zlabel('Pa');
xlabel('X(m)');
ylabel('Y(m)');
title(titbuf);
lim = axis;
M=moviein(nframes);

for n=1:nframes,

  % pressure calculation
  P1 =j.*k.*rc.*Q.*exp(-j.*k.*R1).*exp(j*2*pi*f*n*T)./(4.*pi.*R1);
  P2 =-j.*k.*rc.*Q.*exp(-j.*k.*R2).*exp(j*2*pi*f*n*T)./(4.*pi.*R2);
  P = P1 + P2;
  surfl(X,Y,imag(P));
  shading interp;
  zlabel('Pa');
  xlabel('X(m)');
  ylabel('Y(m)');
  title(titbuf);
  axis(lim);
  M(:,n)=getframe;

end
movie(M,300);
	
