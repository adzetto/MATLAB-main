% Kalman filter - adapted for prognostics example 5/8/99

% sigmaw is measurement std (position only)
% sigmav is process noise (assume white jerk - piecewise const accel)

sigmaw=10;	% temp resolution
sigmav=0.025/3600;		% process noise for Kalman that works well is .05
dt=60;			% 1 hour update in minutes (approx)

% smoother for stable point

R = sigmaw^2;							% 1x1 measurement covariance

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

npts=200;

xku=zeros([3,npts]);	% updated and predicted states
xkp=zeros([3,npts]);	% for plotting

xvp=zeros([1,npts]);	%std's for plotting
xvm=zeros([1,npts]);
vvp=zeros([1,npts]);
vvm=zeros([1,npts]);
avp=zeros([1,npts]);
avm=zeros([1,npts]);

xstdk=zeros([1,npts]);
vstdk=zeros([1,npts]);

Time=zeros([1,npts]);
Peakf=zeros([1,npts]);

%[file path]=uigetfile('*.dat');
%filename=sprintf('%s%s',path,file);
%FP=fopen(filename,'r');
FP=fopen('prognostics.dat','r');
rec=0;

% initialization 
	F = [1 dt .5*dt*dt ; 
	 		0 1  dt;
     		0 0  1]; 		% state transistion matrix
	Q = sigmav^2.*[.25*dt^4 .5*dt^3 .5*dt^2; 
			.5*dt^3  dt^2    dt;
			.5*dt^2  dt      1];	% process covariance


[Time(1) count]=fscanf(FP,'%e',1);
[Peakf(1) count]=fscanf(FP,'%e',1);

Time(1)=Time(1)*60; % time in minutes

xu(1)=Peakf(1);
xku(:,1)=xu;
xkp(:,1)=xu;
Pu=Q;
rec=1;

for k=2:npts,
	[data count]=fscanf(FP,'%e',2);
	if count ~= 2 
		break;
	end;
	rec=rec+1;
	Time(k)=data(1)*60; % T in minutes
	Peakf(k)=data(2);

	dt=(Time(k)-Time(k-1));	% T in minutes
	%dt=(Time(k)-Time(k-1));		% T in seconds
	F = [1 dt .5*dt*dt ; 
 		0 1  dt;
     		0 0  .5]; 		% state transistion matrix
	Q = sigmav^2.*[.25*dt^4 .5*dt^3 .5*dt^2; 
			.5*dt^3  dt^2    dt;
			.5*dt^2  dt      1];	% process covariance

	xp=F*xu;					% predicted
	Pp=F*Pu*F' + Q;
	zp=H*xp;

	ep=Peakf(k) - zp;		% error
	S=H*Pp*H' + R;
	K=Pp*H'*S^(-1);

	xu=xp + K*ep;			% updated
%	Pu=Pp - K*S*K';
 	Pu=(I-K*H)*Pp*(I-K*H)' + K*R*K';	% Joseph form (cleaner)

	xku(:,k)=xu;			% save for plotting
	xkp(:,k)=xp;
	
	xstd=sqrt(Pu(1,1));
	vstd=sqrt(Pu(2,2));
	astd=sqrt(Pu(3,3));

	xvp(k)=xp(1)+xstd;	
	xvm(k)=xp(1)-xstd;	
	vvp(k)=xp(2)+vstd;	
	vvm(k)=xp(2)-vstd;	
	avp(k)=xp(3)+astd;	
	avm(k)=xp(3)-astd;	

	xstdk(k)=xstd;	% save stds for tf variance calc
	vstdk(k)=vstd;
end;

nstart=10; % let filter stabilize

figure(1);
subplot(3,1,1),plot(Time(nstart:rec)./60,xkp(1,nstart:rec),'k',Time(nstart:4:rec)./60,Peakf(nstart:4:rec),'k.',Time(nstart:rec)./60,xvp(nstart:rec),'k:',Time(nstart:rec)./60,xvm(nstart:rec),'k:');
ylabel('Deg C');
subplot(3,1,2),plot(Time(nstart:rec)./60,xkp(2,nstart:rec).*60,'k');
ylabel('Deg C/hr');
subplot(3,1,3),plot(Time(nstart:rec)./60,xkp(3,nstart:rec).*3600,'k');
ylabel('Deg C/hr^2');
xlabel('Time Hours');

% time until failure

xthreshold = 300;	% 150 C for repair, 300 C for failure
tactual=180;

nstart=1;

xkp(2,1)=1;
xkp(2,2)=1;
tf=(xthreshold - xkp(1,:))./xkp(2,:);	

for nn=1:rec-nstart+1,
   if tf(nn) < -6000
      tf(nn)= 20000;
   end;
end;

% we know that the pdf for tf has "std" 1.836 times the Gaussian of (xf-x)/v
% for each tf get the variance of xf

tfp=tf;
tfm=tf;
tstd=zeros(size(tf));
tr=zeros(size(tf));		% survivor
th=zeros(size(tf));		% hazard

for k=nstart:rec,
   dt=tf(k-nstart+1);	% T in minutes
   dtp=dt;
   if dt < 0
      dt=0;					% if you're beyond threshold its clamped
   end;
   
	F = [1 dt .5*dt*dt ; 
 		0 1  dt;
     		0 0  .5]; 		% state transistion matrix
	Q = sigmav^2.*[.25*dt^4 .5*dt^3 .5*dt^2; 
			.5*dt^3  dt^2    dt;
			.5*dt^2  dt      1];	% process covariance
	Pp=F*Pu*F' + Q;				% use same Pu ...

	xfstd=sqrt(Pp(1,1));
	a=xfstd/(vstdk(k)+1e-9);
   tstd(k-nstart+1)=1.836*a;
   
   tr(k-nstart+1)=.5-(1/pi)*atan(-dtp/a);	% survivor
   pdf=(a/pi)/(a^2 + dt^2);
   th(k-nstart+1)=pdf/tr(k-nstart+1);		% hazard
end;
tfp=tf+tstd;
tfm=tf-tstd;
hazmax=max(th);
th=th./hazmax;		% normalize hazard to 1 max

nstable=15;

figure(2);

plot(Time./60,tf./60,'k');	% units of hours, vel & Tf in min
hold on;
plot(Time./60,((tactual*60).*ones(size(Time))-Time)./60,'k--');
plot(Time./60,tfp./60,'k:');
plot(Time./60,tfm./60,'k:');
hold off;
axis([0 200 -50 300]);
ylabel('RUL Hours');
xlabel('Time Hours');

figure(3);
plot(Time(nstart+nstable:rec)./60,tr(1+nstable:rec-nstart+1),'k:');
hold on;
plot(Time(nstart+nstable:rec)./60,th(1+nstable:rec-nstart+1),'k');
hold off;
xlabel('Time Hours');
axis([Time(nstart+nstable)/60 Time(rec)/60 0 1]);
legend('Survivor Probability','Normalized Hazard Rate');

figure(4);

Nbar=20;
alpha=(Nbar-1)/Nbar;
beta=1-alpha;

stability = -xkp(3,:).*xkp(2,:);
for n=2:rec,
   stability(n)=alpha*stability(n-1)+beta*stability(n);	%smoothing
end;

plot(Time(10:rec)./60,stability(10:rec),'k');
ylabel('Stability Metric');
xlabel('Time Hours');

figure(5);
tfc=tf;			% plot clipped tf in 3d - overlap pdf
tmax=300;
tmax=tmax*60;
for n=1:rec-nstart+1,
   if tfc(n) > tmax
      tfc(n)=tmax;
   end;
end;
zz=zeros(size(tf));
plot3(Time(nstart+nstable:rec)./60,tfc(1+nstable:rec-nstart+1)./60,zz(1+nstable:rec-nstart+1),'k');

hold on;

plot3(Time(nstart+nstable:rec)./60,((tactual*60).*ones(size(Time(nstart+nstable:rec)))-Time(nstart+nstable:rec))./60,zz(1+nstable:rec-nstart+1),'k--');

% overlay last pdf
nt=rec;	% last record at failure
a=tstd(nt-nstart+1)/1.836;
t=linspace(-30,tmax/60,200);		% hours
t=t-(Time(rec)/60-tactual);
tm=t.*60;			% minutes
tpdf=(a/pi)./(a^2+(tm-tf(nt-nstart+1)).^2);
tpdfmax1=max(tpdf);
plot3((Time(rec)/60).*ones(size(t)),t,tpdf,'k');

% what is the index into t corresponding to Time(rec)-Time(nt)? (0 hour)
trul=(tactual*60-Time(nt)-t(1)*60)/60;	% hours
ntstart=round(trul/(t(2)-t(1)));
ntstop=round((tmax/60)/(t(2)-t(1)));
for n=ntstart:5:ntstop,
   line([Time(nt)/60 Time(nt)/60],[t(n) t(n)],[0 tpdf(n)],'Color','k');
end;

nt=rec-20;
a=tstd(nt-nstart+1)/1.836;
tpdf=(a/pi)./(a^2+(tm-tf(nt-nstart+1)).^2);
tpdfmax2=max(tpdf);
plot3((Time(nt)/60).*ones(size(t)),t,tpdf,'k');

% what is the index into t corresponding to Time(rec)-Time(nt)? (tinterest hour)
trul=(tactual*60-Time(nt)-t(1)*60)/60;	% hours
ntstart=round(trul/(t(2)-t(1)));
ntstop=round((tmax/60)/(t(2)-t(1)));
for n=ntstart:5:ntstop,
   line([Time(nt)/60 Time(nt)/60],[t(n) t(n)],[0 tpdf(n)],'Color','k');
end;



nt=rec-50;
a=tstd(nt-nstart+1)/1.836;
tpdf=(a/pi)./(a^2+(tm-tf(nt-nstart+1)).^2);
tpdfmax3=max(tpdf);
plot3((Time(nt)/60).*ones(size(t)),t,tpdf,'k');

% what is the index into t corresponding to Time(rec)-Time(nt)? (tinterest hour)
trul=(tactual*60-Time(nt)-t(1)*60)/60;	% hours
ntstart=round(trul/(t(2)-t(1)));
ntstop=round((tmax/60)/(t(2)-t(1)));
for n=ntstart:5:ntstop,
   line([Time(nt)/60 Time(nt)/60],[t(n) t(n)],[0 tpdf(n)],'Color','k');
end;

hold off;
tpdfmax=max([tpdfmax1 tpdfmax2 tpdfmax3]);

axis([Time(1)/60 Time(rec)/60 -50 tmax/60 0 tpdfmax]);
grid on;
view(-30,30);
xlabel('Time Hours');
ylabel('RUL Hours');
zlabel('Probability Density');
title('PDF and Survivor Curves');
