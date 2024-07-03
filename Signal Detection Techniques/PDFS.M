% gauss, chi-1, exponential (chi-2)and Rayleigh pdfs

sigx=1;			% underlying Gaussian std dev
beta=2*sigx^2;
Mx=0;			% x zero mean

npts=20;
ymax=5*sigx;
dy=ymax/npts;
y=dy:dy:ymax;
yg=-ymax+dy:dy:ymax;

pyg=(1./(sigx*sqrt(2*pi))).*exp(-(yg.^2)./beta);
p1y=(1./(sigx.*sqrt(2.*pi.*y))).*exp(-(y./beta));
p2y=(1/beta).*exp(-(y./beta));
pry=(2.*y./beta).*exp(-(y.^2)./beta);

P1y=0;
P2y=0;
Pry=0;
for n=1:npts,
	P1y=P1y+p1y(n)*dy;
	P2y=P2y+p2y(n)*dy;
	Pry=Pry+pry(n)*dy;
end;
Pyg=0;
for n=1:2*npts,
	Pyg=Pyg+pyg(n)*dy;
end;
[P1y P2y Pry Pyg]

figure(1);
plot(yg,pyg,'k');
ylabel('p(y)');
xlabel('y');
hold on;
plot(y,p1y,'k-+');
plot(y,p2y,'k-O');
plot(y,pry,'k-*');
legend('Gaussian','Chi-Sq 1 DOF','Exponential','Rayleigh');
hold off;
axis([-5 5 0 1]);
