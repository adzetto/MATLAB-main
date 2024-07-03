% 2-D FFT for beam response of 8 element line array
ntheta=400;
theta=linspace(-pi+1e-11,pi+1e-11-pi/ntheta,ntheta);
thetadeg=theta*180/pi;

% user input
steerdeg=30;
xmax=32;
ymax=32;		% 32m by 32m 2-d fft
npts=256;
nelements=4;	% works for 2,4,8,16 elements
L=2;			   % 2 meters
f=600;			% 300 Hz	

buff=sprintf('%dHz %d Sensors L=%d Steer=%d deg',f,nelements,L,steerdeg);


steer=steerdeg*pi/180;
j=sqrt(-1);

c=344;		% wave speed
lambda=c/f;
%L=lambda/2;
k=2*pi/lambda;

% theoretical directivity including steering
D=L.*sin((k*L/2).*(sin(theta)-sin(steer)))./((k*L/2).*(sin(theta)-sin(steer)));

figure(1);
plot(thetadeg,abs(D),'k');
xlabel('\theta');
ylabel('Beam Response');
axis([-180 180 0 2.5]);

% Get directivity index including direction
intD=0;
for n=1:ntheta,
	intD=intD+D(n)^2;
end;
intD=intD/ntheta;

nsteer=round(steer*ntheta/(2*pi))+ntheta/2;
DIest=(D(nsteer)^2)/intD	% Directivity index (linear)
DI=2*L/lambda	% Burdic's approximation

% 2-D FFT beamforming representation

dx=xmax/npts;
dy=ymax/npts;
L2=L/2;
nl=round(L2/dx);	% length of line array
xsp=L/nelements     % spacing elements
nsp=round(xsp/dx);	% pixel spacing of line

% load the line source on the grid
S=zeros([npts npts]);
for i=-nl:nsp:nl-1,
	S(128,128+i)=exp(j*k*i*dx*sin(steer));
end;
figure(3);
image(64.*abs(S));
colormap(gray);
xlabel('Input X Bin');
ylabel('Input Y Bin');
title(buff);

Sk=fftshift(fft2(S));

Skl=abs(Sk);
maxS=max(max(Skl));		% scaling for figure only

figure(4);
imagesc(Skl);
colormap(gray);
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

Sscale=dx*npts*L/(xmax*nelements); 
% L = dx*npts*L/xmax
% actual wavenumber spectrum scaling 
% kmax = 2*pi/xmax, dk=kmax/npts
Dkk=zeros(size(theta));
for n=1:ntheta,
	Dkk(n)=Sscale*Sk(nky(n),nkx(n));
end;

figure(5);
plot(thetadeg,abs(D),'k',thetadeg,abs(Dkk),'k.');
title(buff);
xlabel('Degrees');
ylabel('Beam Response');
legend('Theory','2-D K-Space Resp');

figure(2);
polar(theta,abs(Dkk),'k');
title(buff);
view(90,-90);

