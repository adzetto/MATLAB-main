% point source field planes
% reconstruction of quadrapole field at z1
% from measurement at z2
% transfer fcn method

f=700;
%f=1024;     % increase f restores resolution
c=344;
k=2*pi*f/c;

z1=.001;	% reconstruction planes
z2=5;       % measurement plane
%z2=1.0;

xwidth=10;
ywidth=10;
npts=64;
npts=32;
dx=xwidth/npts;
dy=ywidth/npts;

% plane grid
[X Y]=meshgrid(-xwidth/2:dx:xwidth/2-dx,-ywidth/2:dy:ywidth/2-dy);

R1=(z1^2 + X.^2 + Y.^2).^.5;	% ranges to x-y grid at z1
R2=(z2^2 + X.^2 + Y.^2).^.5;

j=sqrt(-1);
P1=exp(-j.*k.*R1)./(4.*pi.*R1);	% complex field on plane 1
P2=exp(-j.*k.*R2)./(4.*pi.*R2);


% wavenumber params
kxs=npts/xwidth;
kys=npts/ywidth;
dkx=2*pi*kxs/npts;
dky=2*pi*kys/npts;

[KX KY]=meshgrid(-pi*kxs:dkx:pi*kxs-dkx,-pi*kys:dky:pi*kys-dky);

P1F=fftshift(fft2(P1));
P2F=fftshift(fft2(P2));

% point source transfer function
HF=P1F./P2F;	

% quadrapole sources
xqp=[1.2 -1 -1.8 1.5];
yqp=[1 1.3 -1.5 -.5];

%xqp=[1 -1 -1 1];
%yqp=[1 1 -1 -1];   % equally spaces sources

pzq=[pi 0 pi 0];

QP1=zeros(size(P1));	% actual field
QP2=zeros(size(P2));	% measurement field
RQP1=zeros(size(P1));	% reconstructed field

for n=1:4,
	R1=(z1^2 + (X-xqp(n)).^2 + (Y-yqp(n)).^2).^.5;
	R2=(z2^2 + (X-xqp(n)).^2 + (Y-yqp(n)).^2).^.5;

	QP1=QP1+exp(j*pzq(n)).*exp(-j.*k.*R1)./(4.*pi.*R1);
	QP2=QP2+exp(j*pzq(n)).*exp(-j.*k.*R2)./(4.*pi.*R2);
end;

QP2F=fftshift(fft2(QP2));	% measuremnt wavenumber field

QP1F=HF.*QP2F;	% reconstruction
RQP1=ifft2(fftshift(QP1F));


% check
figure(1);

set(gcf,'color','w');

title1=sprintf('Real Pressure at %.2f m',z2);
subplot(2,2,3),pcolor(X,Y,real(QP2));
%axis([-6 6 -6 6]);
set(gca,'color','w');
set(gca,'xcolor','k');
set(gca,'ycolor','k');
set(gca,'fontsize',8);
colormap(gray);
shading interp;
title(title1,'color','k');
xlabel('x[m]','color','k');
ylabel('y[m]','color','k');

title2=sprintf('Reconstruction at %.2f m',z1);
subplot(2,2,1),pcolor(X,Y,real(RQP1));
%axis([-6 6 -6 6]);
set(gca,'color','w');
set(gca,'xcolor','k');
set(gca,'ycolor','k');
set(gca,'fontsize',8);
shading interp;
title(title2,'color','k');
xlabel('x[m]','color','k');
ylabel('y[m]','color','k');

title3=sprintf('Wavenumber Spectrum at %.2f m',z2);
subplot(2,2,4),pcolor(KX,KY,abs(QP2F));
%axis([-25 25 -25 25]);
set(gca,'color','w');
set(gca,'xcolor','k');
set(gca,'ycolor','k');
set(gca,'fontsize',8);
colormap(gray);
shading interp;
title(title3,'color','k');
xlabel('Kx','color','k');
ylabel('Ky','color','k');

title4=sprintf('Wavenumber Spectrum at %.2f m',z1);
subplot(2,2,2),pcolor(KX,KY,abs(QP1F));
%axis([-25 25 -25 25]);
set(gca,'color','w');
set(gca,'xcolor','k');
set(gca,'ycolor','k');
set(gca,'fontsize',8);
colormap(gray);
shading interp;
title(title4,'color','k');
xlabel('Kx','color','k');
ylabel('Ky','color','k');

