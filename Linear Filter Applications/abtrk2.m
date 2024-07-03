npts=100;
xtrue=zeros(size([1,npts]));
vtrue=zeros(size([1,npts]));
atrue=zeros(size([1:npts]));
xtrue(1)=0;
vtrue(1)=0;
atrue(1)=0;

agravity=9.801;
drag=.25;
arocket=15;
burnout=50;
dt=.1;
t=dt:dt:npts*dt;

for k=2:burnout,
	atrue(k)=arocket-agravity-drag*vtrue(k-1);
	vtrue(k)=vtrue(k-1)+atrue(k-1).*dt;
	xtrue(k)=xtrue(k-1)+vtrue(k-1).*dt+.5.*atrue(k-1).*dt^2;
end
for k=burnout:100,
	atrue(k)=-agravity-drag*vtrue(k-1);
	vtrue(k)=vtrue(k-1)+atrue(k-1).*dt;
	xtrue(k)=xtrue(k-1)+vtrue(k-1).*dt+.5.*atrue(k-1).*dt^2;
end

xmeas=zeros(size([1,npts]));
vmeas=zeros(size([1,npts]));
ameas=zeros(size([1:npts]));

sigmaw = 3;

xmeas(1)=xtrue(1)+randn.*sigmaw;
vmeas(1)=xmeas(1);
ameas(1)=vmeas(1);
burston = 0;

for k=2:npts,
	if k>70 & k<80 & burston == 1
		sigmaw=40;
	else
		sigmaw=3;
	end
	xmeas(k)=xtrue(k)+randn.*sigmaw;
	vmeas(k)=(xmeas(k)-xmeas(k-1))./dt;
	ameas(k)=(vmeas(k)-vmeas(k-1))./dt;
end

xpred=zeros(size([1,npts]));
vpred=zeros(size([1,npts]));
apred=zeros(size([1:npts]));

xupdat=zeros(size([1,npts]));
vupdat=zeros(size([1,npts]));
aupdat=zeros(size([1:npts]));

sigmav=100;

xupdat(1)=xmeas(1);
vupdat(1)=0;
aupdat(1)=0;
trkburst = 0;

for k=2:npts,

	if k>70 & k<80 & trkburst == 1
		sigmaw=40;
	else
		sigmaw=3;
	end
	l=sigmav*dt^2/sigmaw;
	b=.25*(l^2+4*l-l*sqrt(l^2+8*l));
	a=sqrt(2*b)-b/2;
	g=b^2/a;

	xpred(k)=xupdat(k-1) + dt*vupdat(k-1) + (.5*dt^2)*aupdat(k-1);
	vpred(k)=vupdat(k-1) + dt*aupdat(k-1);

%	xpred(k)=xupdat(k-1) + dt.*vupdat(k-1);  % a-b only
%	vpred(k)=vupdat(k-1);

	apred(k)=aupdat(k-1);

	error=xmeas(k)-xpred(k);

	xupdat(k)=xpred(k)+a*error;
	vupdat(k)=vpred(k)+b*error/dt;
	aupdat(k)=apred(k)+g*error/(2*dt^2); 
end

figure(1);
subplot(3,1,1),plot(t,xmeas,'ok',t,xpred,'k');
subplot(3,1,2),plot(t,vpred,'k');
subplot(3,1,3),plot(t,apred,'k');



