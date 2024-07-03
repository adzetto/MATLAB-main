% 2-D FFT for beam response of 8 element line array
ntheta=400;
theta=linspace(-pi+1e-11,pi+1e-11-pi/ntheta,ntheta);
thetadeg=theta*180/pi;

% user input
steerdeg=30;
xmax=32;
ymax=32;		% 32m by 32m 2-d fft
npts=256;
nxels=8;		% works for 2,4,8,16 elements
nyels=8;
Lx=2;			% 2 meters
Ly=2;
f=300;			% 300 Hz	

buff=sprintf('%dHz %d Sensors %dx%dm Steer=%d deg',f,nxels*nyels,Lx,Ly,steerdeg);

steer=steerdeg*pi/180;
j=sqrt(-1);

c=344;		% wave speed
lambda=c/f;
k=2*pi/lambda;

% theoretical directivity including steering
Dx=Lx.*sin((k*Lx/2).*(sin(theta)-sin(steer)))./((k*Lx/2).*(sin(theta)-sin(steer)));
Dy=Ly.*sin((k*Ly/2).*(cos(theta)-cos(steer)))./((k*Ly/2).*(cos(theta)-cos(steer)));
D=Dx.*Dy;


% Get directivity index including direction
intD=0;
for n=1:ntheta,
	intD=intD+D(n)^2;
end;
intD=intD/ntheta;

nsteer=round(steer*ntheta/(2*pi))+ntheta/2;
DIest=(D(nsteer)^2)/intD	% Directivity index (linear)

% 2-D FFT beamforming representation

dx=xmax/npts;
dy=ymax/npts;
Lx2=Lx/2;
Ly2=Ly/2;
nlx=round(Lx2/dx);	% length of x line array in pixels
xsp=Lx/nxels;	   	% x spacing elements
nspx=round(xsp/dx);	% pixel spacing of x line
nly=round(Ly2/dy);	% length of y line array in pixels
ysp=Ly/nyels;	   	% y spacing elements
nspy=round(ysp/dy);	% pixel spacing of x line

% load the line source on the grid
S=zeros([npts npts]);
for ix=-nlx:nspx:nlx-1,
  for iy=-nly:nspy:nly-1,
	S(128+iy,128+ix)=exp(j*k*ix*dx*sin(steer))*exp(j*k*iy*dy*cos(steer));
  end;
end;
figure(2);
image(64.*abs(S));
colormap(gray);
xlabel('Input X Bin');
ylabel('Input Y Bin');
title(buff);

Sk=fftshift(fft2(S));

Skl=abs(Sk);
maxS=max(max(Skl));		% scaling for figure only

figure(3);
image((256/maxS).*Skl);
colormap(gray);
%colorbar;
xlabel('FFT Kx Bin');
ylabel('FFT Ky Bin');
title(buff);

% define circle for k
dkx=2*pi/xmax;
dky=2*pi/ymax;

nkx=npts/2 + round((k/dkx).*sin(theta)+1);	% shift 1 pixel
nky=npts/2 + round((k/dky).*cos(theta)+1);

hold on;
plot(nkx,nky,'r--');
hold off;

Sscale=(dx*npts*Lx/(xmax*nxels))*(dy*npts*Ly/(ymax*nyels));
% L = dx*npts*L/xmax
% actual wavenumber spectrum scaling 
% kmax = 2*pi/xmax, dk=kmax/npts
Dkk=zeros(size(theta));
for n=1:ntheta,
	Dkk(n)=Sscale*Sk(nky(n),nkx(n));
end;

figure(4);
plot(thetadeg,abs(D),'k',thetadeg,abs(Dkk),'k.');
title(buff);
xlabel('Degrees');
ylabel('Beam Response');
legend('Theory','2-D K-Space Resp');

figure(1);
polar(theta,abs(Dkk),'k');
title(buff);
view(90,-90);


% Get directivity index including direction
intD=0;
for n=1:ntheta,
	intD=intD+abs(Dkk(n)^2);
end;
intD=intD/ntheta;

DIestkk=abs(Dkk(nsteer)^2)/intD	% Directivity index (linear)

