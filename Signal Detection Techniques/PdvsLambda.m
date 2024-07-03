% Pd vs detection for time-averaged RMS detector

fs=4096;        % sample rate
Tmax=1;        % max recording time in sec
TN1=.01;          % sig time average in sec
TN2=.1;         % bkg time average in sec

% generate unity variance white noise signal
npts=Tmax*fs;
T=[1:npts]./fs;
nz=randn(1,npts);

% lets put a simple pulse in there 200 msec
Tpms=200;
Poffset=1000;
amp=1;
Plen=round((Tpms/1000)*fs); % duration in samples
nz(Poffset:Poffset+Plen)=nz(Poffset:Poffset+Plen)+amp.*ones([1,Plen+1]);

nzsq=nz.^2;     % noise squared

a1=(TN1*fs-1)/(TN1*fs); % exp integrators
b1=1-a1;
a2=(TN2*fs-1)/(TN2*fs);
b2=1-a2;

sig=zeros(size(nz));    % signal is short time average
bkg=zeros(size(nz));    % background is long time average

sigsq=zeros(size(nz));  % for mean square calc
bkgsq=zeros(size(nz));

% calculate mean square with exponential integration window
sigsq(1)=1; % start them out a the stdev
bkgsq(1)=1;
for n=2:npts,
    sigsq(n)=a1*sigsq(n-1)+b1*nzsq(n);
    bkgsq(n)=a2*bkgsq(n-1)+b2*nzsq(n);
end;

sig=sqrt(sigsq);    % get back the RMS from mean-square
bkg=sqrt(bkgsq);

figure(1);
pltpts=4096;
subplot(3,1,1),plot(T(1:pltpts),nz(1:pltpts),'k');
text(.6,amp,sprintf('\\sigma = 1.00  A=%.1f',amp));
ylabel('Input');
title('RMS Averaging of Signals in ZMG Noise');
%axis([0, T(pltpts),-5, amp+2]);

subplot(3,1,2),plot(T(1:pltpts),sig(1:pltpts),'k');
text(.6,.8*amp,sprintf('\\tau = %.3f sec',TN1));
ylabel('Fast RMS');
%axis([0, T(pltpts),0, amp+2]);

subplot(3,1,3),plot(T(1:pltpts),bkg(1:pltpts),'k');
text(.6,.8*amp,sprintf('\\tau = %.3f sec',TN2));
ylabel('Slow RMS');
%axis([0, T(pltpts),0, amp+2]);
xlabel('Sec');

% Ok now lets "brute force" the Pd by counting samples above the threshold
% and sweep the threshold up from zero on the amplitude data.

nzamp=sqrt(nzsq);   % instantaneous amplitude

ntpts=100;
gammax=amp+2;
gamma=[ntpts:-1:1].*gammax./ntpts;

detcntamp=zeros(size(gamma));
detcntsig=zeros(size(gamma));
detcntbkg=zeros(size(gamma));

for k=1:ntpts,
    for n=1:npts,
        if nzamp(n)>gamma(k)
            detcntamp(k)=detcntamp(k)+1;
        end;
        if sig(n)>gamma(k)
            detcntsig(k)=detcntsig(k)+1;
        end;
        if bkg(n)>gamma(k)
            detcntbkg(k)=detcntbkg(k)+1;
        end;
    end;
end;
% normalize to get Pd
detcntamp=detcntamp./npts;
detcntsig=detcntsig./npts;
detcntbkg=detcntbkg./npts;

figure(2);
plot(gamma,detcntamp,'k',gamma,detcntsig,'k-.',gamma,detcntbkg,'k:');
legend('Input Amp','Fast RMS','Slow RMS','Location','NorthWest');
xlabel('Absolute Detection Threshold {\Lambda}'); % cap first letter for uppercase Greek
ylabel('Estimated Pd');
% show the plot with the threshold axis going backwards (rotate 180 and
% flip 180)
view(180,-90);          % standard plot view is view(0,90)
title(sprintf('Noise \\sigma_N=1.00 Signal Amp = %.2f',amp));
