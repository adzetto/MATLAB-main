% pulse generator and integrals

% T/(2W) gives number of harmonics in freq domain!
T=8500;	% sample # pulse repeats
W=4;		% pulse width in samples
N=(T/W);		% T/N is pulse width in samples
tau=200.5;	% pulse delay
% if T/W is odd -> even harmonics
% if T/W is even -> odd harmonics


npts=400;
n=1:npts;

y=(1/N).*sin(pi.*N.*(n-tau)./T)./(sin(pi.*(n-tau)./T));

for nn=1:npts,
   if y(nn)>.97
      y(nn)=1;
   end;
   if y(nn)<-.97
      y(nn)=-1;
   end;
end;

pulse=y;

% integrate the waveforms
step=zeros(size(y));
ramp=zeros(size(y));
quad=zeros(size(y));
step(1)=pulse(1);
ramp(1)=step(1);
quad(1)=ramp(1);
for nn=2:npts,
   step(nn)=step(nn-1)+pulse(nn);
   ramp(nn)=ramp(nn-1)+step(nn);
   quad(nn)=quad(nn-1)+ramp(nn);
end;

% differentiate the waveforms
pulsed=diff(pulse);
stepd=diff(step);
rampd=diff(ramp);
quadd=diff(quad);


figure(1);
subplot(4,2,1),plot(n,pulse,'k');
title('F[t]');
subplot(4,2,3),plot(n,step,'k');
subplot(4,2,5),plot(n,ramp,'k');
subplot(4,2,7),plot(n,quad,'k');
xlabel('t');

subplot(4,2,2),plot(n(1:npts-1),pulsed,'k');
title('dF[t]/dt');
subplot(4,2,4),plot(n(1:npts-1),stepd,'k');
axis([0 400 -1 1]);
subplot(4,2,6),plot(n(1:npts-1),rampd,'k');
subplot(4,2,8),plot(n(1:npts-1),quadd,'k');
xlabel('t');