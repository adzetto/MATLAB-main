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

% rocket Newtonian motion

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

%s0=3;
s0=0.1;         % low measurement noise
sigmaw = s0;	% measurement noise variance

xmeas(1)=xtrue(1)+randn.*sigmaw;
vmeas(1)=xmeas(1);
ameas(1)=vmeas(1);
burston = 0;		% flag for noise burst

for k=2:npts,
	if k>70 && k<80 && burston == 1
		sigmaw=40;
	else
		sigmaw=s0;
	end
	xmeas(k)=xtrue(k)+randn.*sigmaw;
	vmeas(k)=(xmeas(k)-xmeas(k-1))./dt;	% these are for comparison only
	ameas(k)=(vmeas(k)-vmeas(k-1))./dt;
end

xpred=zeros(size([1,npts]));
vpred=zeros(size([1,npts]));
apred=zeros(size([1:npts]));

xupdat=zeros(size([1,npts]));
vupdat=zeros(size([1,npts]));
aupdat=zeros(size([1:npts]));

%sigmav=13;		% high process noise std dev (13 optimal for abg)
sigmav=1;		% low process noise std dev (13 optimal for abg)

% alpha-beta-gamma tracker

xupdat(1)=xmeas(1);
vupdat(1)=0;
aupdat(1)=0;
trkburst = 0;

for k=2:npts,

	if k>70 && k<80 && trkburst == 1
		sigmaw=40;
	else
		sigmaw=s0;
	end
	l=sigmav*dt^2/sigmaw;				% maneuverability
	b=.25*(l^2+4*l-l*sqrt(l^2+8*l));	% beta
	a=sqrt(2*b)-b/2;					% alpha
	g=b^2/a;							% gamma

	xpred(k)=xupdat(k-1) + dt.*vupdat(k-1) + (.5*dt^2).*aupdat(k-1);
	vpred(k)=vupdat(k-1) + dt.*aupdat(k-1);
	apred(k)=aupdat(k-1);

%	xpred(k)=xupdat(k-1) + dt.*vupdat(k-1);  % a-b only
%	vpred(k)=vupdat(k-1);

	error=xmeas(k)-xpred(k);

	xupdat(k)=xpred(k)+a*error;
	vupdat(k)=vpred(k)+b*error/dt;
	aupdat(k)=apred(k)+g*error/(2*dt^2); 
end

figure(1);
subplot(3,1,1), plot(t,xtrue,'k',t,xmeas,'.k',t,xpred,'--k');
ylabel('Position m');
xlabel('Seconds');
legend('True Position','Measurements','Tracked','Location','NorthWest');
title('\alpha-\beta-\gamma Tracker');

subplot(3,1,2), plot(t,vtrue,'k',t,vmeas,'.k',t,vpred,'--k');
ylabel('Velocity m/s');
xlabel('Seconds');
legend('True Velocity','Measurements','Tracked','Location','SouthWest');

subplot(3,1,3), plot(t,atrue,'k',t,apred,'--k');
ylabel('Acceleration m/s2');
xlabel('Seconds');
legend('True Acceleration','Tracked','Location','SouthWest');

% Kalman filter
% sigmaw is measurement std (position only)
% sigmav is process noise (assume white jerk - piecewise const accel)

sigmav=1;		% process noise for Kalman that works well
%sigmav=13;		% process noise for Kalman that works well


F = [1 dt .5*dt*dt ; 
	 0 1  dt;
     0 0  1]; 							% state transistion matrix

R = sigmaw^2;							% 1x1 measurement covariance

Q = sigmav^2.*[.25*dt^4 .5*dt^3 .5*dt^2; 
				.5*dt^3  dt^2    dt;
				.5*dt^2  dt      1];	% process covariance

H = [1 0 0];							% measurement matrix

xu=zeros([3 1]);	% a posteriori updated state vector
xp=zeros([3 1]);	% a priori predicted state vector
Pu=zeros([3 3]);	% updated state error covariance
Pp=zeros([3 3]);	% predicted state covariance
K=zeros([3 1]);		% Kalman gain vector
S=zeros([1 1]);
I=zeros([3 3]);
for k=1:3,
	I(k,k)=1;
end;

xku=zeros([3,npts]);	% updated and predicted states
xkp=zeros([3,npts]);	% for plotting
zp=zeros([1,npts]);		% predicted measurements
ep=zeros([1,npts]);		% measurement prediction error

% initialization 
xu(1)=xmeas(1);
xku(:,1)=xu;
xkp(:,1)=xu;
Pu=Q;

for k=2:npts,

	if k>70 && k<80 && trkburst == 1
		sigmaw=40;
	else
		sigmaw=s0;
	end
	R = sigmaw^2;		% 1x1 measurement covariance


	xp=F*xu;
	Pp=F*Pu*F' + Q;
	zp=H*xp;
	ep=xmeas(k) - zp;
	S=H*Pp*H' + R;
	K=Pp*H'*S^(-1);
	xu=xp + K*ep;
%	Pu=Pp - K*S*K';
    Pu=(I-K*H)*Pp*(I-K*H)' + K*R*K';	% Joseph form (cleaner)
    
	xku(:,k)=xu;
	xkp(:,k)=xp;

end;

figure(2);
subplot(3,1,1), plot(t,xtrue,'k',t,xmeas,'.k',t,xkp(1,:),'--k');
ylabel('Position m');
xlabel('Seconds');
legend('True Position','Measurements','Tracked','Location','NorthWest');
title('Kalman filter tracker');

subplot(3,1,2), plot(t,vtrue,'k',t,vmeas,'.k',t,xkp(2,:),'--k');
ylabel('Velocity m/s');
xlabel('Seconds');
legend('True Velocity','Measurements','Tracked','Location','SouthWest');

subplot(3,1,3), plot(t,atrue,'k',t,xkp(3,:),'--k');
ylabel('Acceleration m/s2');
xlabel('Seconds');
legend('True Acceleration','Tracked','Location','SouthWest');

