% 2-D FFT for beam response of line array
ntheta=400;
theta=linspace(-pi+1e-11,pi+1e-11-pi/ntheta,ntheta);
thetadeg=theta*180/pi;

steerdeg=30;
steer=steerdeg*pi/180;
j=sqrt(-1);

L=2;		% 2 meters
f=300;		% 300 Hz	
c=344;		% wave speed
lambda=c/f;
k=2*pi/lambda

% theoretical directivity including steering
D=L.*sin((k*L/2).*(sin(theta)-sin(steer)))./((k*L/2).*(sin(theta)-sin(steer)));

figure(1);
plot(thetadeg,abs(D),'k');
ylabel('Beam Response');
%xlabel('q');  % this is VERY old school...
%set(get(gca,'xlabel'),'FontName','Symbol');
xlabel('\theta');
axis([-180 180 0 2.5]);
figure(2);
polar(theta,abs(D),'k');
view(90,-90);

intD=0;
for n=1:ntheta,
	intD=intD+D(n)^2;
end;
intD=intD/ntheta;

nsteer=round(steer*ntheta/(2*pi))+ntheta/2;
DIest=(D(nsteer)^2)/intD	% Directivity index (linear)
DI=2*L/lambda	% Burdic's approximation

% 2-D FFT beamforming representation
% define grid size:
xmax=32;
ymax=32;	% 32m by 32m
npts=256;
dx=xmax/npts;
dy=ymax/npts;
L2=L/2;
nl=round(L2/dx);	% length of line array

% load the line source on the grid
S=zeros([npts npts]);
for i=-nl:nl-1,
	S(128,128+i)=exp(j*k*i*dx*sin(steer));
end;

Sk=fftshift(fft2(S));

Skl=abs(Sk);
maxS=max(max(Skl));		% scaling for figure only

figure(3);
image(256.*abs(S));
shading interp;
colormap(gray);
xlabel('FFT input x bin');
ylabel('FFT input y bin');

figure(4);
image((256/maxS).*abs(Skl));
colormap(gray);
shading interp;
xlabel('FFT output x bin');
ylabel('FFT output y bin');


% define circle for k
dkx=2*pi/xmax;
dky=2*pi/ymax;

nkx=npts/2 + round((k/dkx).*sin(theta)+1);
nky=npts/2 + round((k/dky).*cos(theta)+1);

hold on;
plot(nkx,nky,'r--');
hold off;

Sscale=xmax/npts;		% actual wavenumber spectrum scaling 
						% kmax = 2*pi/xmax, dk=kmax/npts
Dkk=zeros(size(theta));
for n=1:ntheta,
	Dkk(n)=Sscale*Sk(nky(n),nkx(n));
end;

% This is a check
figure(5);
plot(thetadeg,abs(D),'k',thetadeg,abs(Dkk),'ko');
axis([-180 180 0 2.5]);
xlabel('\theta');
ylabel('FFT Response for K');

