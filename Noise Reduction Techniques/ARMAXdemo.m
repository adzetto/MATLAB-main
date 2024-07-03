% ARMAX demo

fs=1024;
f0=50;

H0=20;
del=1;		% delay must be at least 1

npts=1024;
nwts=10;

BD=zeros([1 nwts]);
BD(1)=1e-1;

n=1:npts;
arg=2*pi*f0/fs;
nt=sin(arg.*n);
yt=zeros(size(nt));
ht=zeros(size(nt));
ut=zeros(size(nt));
utp=zeros([1 npts+del]);
up=zeros([1 del]);

for n=1:nwts
   yt(n)=nt(n);
end;
pow=BD(1)^2;

for n=nwts+1:npts-1,
   Gn=-BD./H0;			% apply inverse gain scale here
  % Gn(1)=(1-BD(1))/H0;
   Gd=BD;
   
   sn=Gn(1)*yt(n);
   sd=0;
   
   % control signal u
   for k=1:nwts,
      sn=sn+Gn(k)*yt(n-k);
      sd=sd+Gd(k)*ut(n-k);
   end;
   ut(n)=sn-sd;
   
   % prediction d-step ahead
   for q=1:del,
      sp=0;
      sn=Gn(q)*yt(n);
      sd=0;
      for k=1:nwts-q,
         sn=sn+Gn(k+q)*yt(n-k);
         sd=sd+Gd(k+q)*ut(n-k);
         if k < q
            sp=sp+Gd(k)*up(q-k);
         end;
         up(q)=sn-sd-sp;
      end;
   end;
   utp(n+del)=up(del);
   
   % whiten control signal
   uhat=0;
   for k=1:nwts,
      uhat=uhat-BD(k)*ut(n-k);
   end;
   
   murel=0.005;     % not too fast or slow!
   pow=.95*pow + .05*std(ut(n-nwts:n))^2;
   mu=1/(pow*nwts);
   err=ut(n)-uhat;
   for k=1:nwts,
      BD(k) = BD(k) - 2*mu*murel*err*ut(n-k);
   end;
 
   % plant output
   ht(n+1)=H0*utp(n-del);		% simulated gain and delay here
   yt(n+1)=nt(n+1)-ht(n+1);
   
end;

figure(1);
plot(yt,'k');
ylabel('n_t');
xlabel('Iteration');
title('ARMAX Active Attenuation Time Response');