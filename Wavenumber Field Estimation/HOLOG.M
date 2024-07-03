% point source field planes

f=700;
c=344;
k=2*pi*f/c;

z1=.0001;	% field planes
z2=1;
z3=50;

xwidth=10;
ywidth=10;
npts=64;
dx=xwidth/npts;
dy=ywidth/npts;

% plane grid
[X Y]=meshgrid(-xwidth/2:dx:xwidth/2-dx,-ywidth/2:dy:ywidth/2-dy);

R1=(z1^2 + X.^2 + Y.^2).^.5;	% ranges to x-y grid at z1
R2=(z2^2 + X.^2 + Y.^2).^.5;
R3=(z3^2 + X.^2 + Y.^2).^.5;

j=sqrt(-1);
P1=exp(-j.*k.*R1)./(4.*pi.*R1);	% complex field on plane 1
P2=exp(-j.*k.*R2)./(4.*pi.*R2);
P3=exp(-j.*k.*R3)./(4.*pi.*R2);


% wavenumber params
kxs=npts/xwidth;
kys=npts/ywidth;
dkx=2*pi*kxs/npts;
dky=2*pi*kys/npts;

[KX KY]=meshgrid(-pi*kxs:dkx:pi*kxs-dkx,-pi*kys:dky:pi*kys-dky);

P1F=fftshift(fft2(P1));
P2F=fftshift(fft2(P2));
P3F=fftshift(fft2(P3));

% check
figure(1);

set(gcf,'color','w');

title1=sprintf('Real Pressure at %.2f m',z1);
subplot(3,2,1),pcolor(X,Y,real(P1));
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

title2=sprintf('|G[kx|ky]| at %.2f m',z1);
subplot(3,2,2),pcolor(KX,KY,abs(P1F));
%axis([-25 25 -25 25]);
set(gca,'color','w');
set(gca,'xcolor','k');
set(gca,'ycolor','k');
set(gca,'fontsize',8);
shading interp;
title(title2,'color','k');
xlabel('Kx','color','k');
ylabel('Ky','color','k');

% data at z2
title1=sprintf('Real Pressure at %.2f m',z2);
subplot(3,2,3),pcolor(X,Y,real(P2));
%axis([-6 6 -6 6]);
set(gca,'color','w');
set(gca,'xcolor','k');
set(gca,'ycolor','k');
set(gca,'fontsize',8);
shading interp;
title(title1,'color','k');
xlabel('x[m]','color','k');
ylabel('y[m]','color','k');

title2=sprintf('|G[kx|ky]| at %.2f m',z2);
subplot(3,2,4),pcolor(KX,KY,abs(P2F));
%axis([-25 25 -25 25]);
set(gca,'color','w');
set(gca,'xcolor','k');
set(gca,'ycolor','k');
set(gca,'fontsize',8);
shading interp;
title(title2,'color','k');
xlabel('Kx','color','k');
ylabel('Ky','color','k');

%data at z3
title1=sprintf('Real Pressure at %.2f m',z3);
subplot(3,2,5),pcolor(X,Y,real(P3));
%axis([-6 6 -6 6]);
set(gca,'color','w');
set(gca,'xcolor','k');
set(gca,'ycolor','k');
set(gca,'fontsize',8);
shading interp;
title(title1,'color','k');
xlabel('x[m]','color','k');
ylabel('y[m]','color','k');

title2=sprintf('|G[kx|ky]| at %.2f m',z3);
subplot(3,2,6),pcolor(KX,KY,abs(P3F));
%axis([-25 25 -25 25]);
set(gca,'color','w');
set(gca,'xcolor','k');
set(gca,'ycolor','k');
set(gca,'fontsize',8);
shading interp;
title(title2,'color','k');
xlabel('Kx','color','k');
ylabel('Ky','color','k');

