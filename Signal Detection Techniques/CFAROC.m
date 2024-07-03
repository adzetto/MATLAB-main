% ROC curves

Nrms=41;    % fast RMS example
%Nrms=410;   % slow RMS example
%Nrms=4;     % raw data

% noise stats
Mn = 1;         % mean level is also std dev for noise
sign = Mn/sqrt(Nrms);    % std dev

% range of interest
dx=0.2;
x=(-5+Mn:dx:12+Mn);

gpdf=(1/(sign*sqrt(2*pi))).*exp((-.5/(sign*sign)).*(x-Mn).^2);

% complimentary erfc
z0snr=(x-Mn)./(sqrt(2)*sign);
Pd0snr=.5*erfc(z0snr);

% verify by integrating Gaussian
[nrow, npts]=size(x);
gsum=zeros(size(x));
for k=1:npts,
    for n=k:npts,
        gsum(k)=gsum(k)+gpdf(n);
    end;
end;
gsum=gsum.*dx;  % dx term for int scaling

% this is just a check for the scaling of the erfc
figure(1);
plot(x,Pd0snr,'k',x,gsum,'k:');
view(180,-90);

% setup signal Pd
Mns1=sqrt(Mn^2+Mn^2);
z1snr=(x-Mns1)./(sqrt(2)*sign);
Mns2=sqrt((2*Mn)^2+Mn^2);
z2snr=(x-Mns2)./(sqrt(2)*sign);
Mns4=sqrt((4*Mn)^2+Mn^2);
z4snr=(x-Mns4)./(sqrt(2)*sign);
Mns10=sqrt((10*Mn)^2+Mn^2);
z10snr=(x-Mns10)./(sqrt(2)*sign);

Pd1snr=.5*erfc(z1snr);
Pd2snr=.5*erfc(z2snr);
Pd4snr=.5*erfc(z4snr);
Pd10snr=.5*erfc(z10snr);

% this is just a check for the scaling of the erfc
figure(2);
semilogy(x,Pd0snr,'k',x,Pd1snr,'k-^',x,Pd2snr,'k-.',x,Pd4snr,'k:',x,Pd10snr,'k--');
legend('SNR 0','SNR 1','SNR 2','SNR 4','SRN 10','Location','SouthEast');
title(sprintf('Pd vs SNR for %d RMS Averages in Window',Nrms));
ylabel('Pd');
xlabel('Absolute Detection threshold \Lambda');
axis([-4 14 1e-8 2]);
view(180,-90);
grid on;

% now we can make a ROC curve

figure(3);
loglog(Pd0snr,Pd1snr,'k-^',Pd0snr,Pd2snr,'k-.',Pd0snr,Pd4snr,'k:',Pd0snr,Pd10snr,'k--');
%view(180,-90);
legend('SNR 1','SNR 2','SNR 4','SNR 10','Location','SouthWest');
title(sprintf('Pd vs Pfa ROC Curves for %d RMS Averages',Nrms));
ylabel('Pd');
xlabel('Pfa');
axis([1e-9 1 1e-2 2]);
grid on;

% lets make a plot of Pfa vs threshold using erfcinv

pfarng=logspace(-12,0,100);
% for the SNR0 case
gamrng=Mn.*ones(size(pfarng)) + (sign*sqrt(2)).*erfcinv(2.*pfarng);

figure(4);
semilogx(pfarng,gamrng,'k');
ylabel('Threshold \Lambda');
xlabel('P_{fa}');
title(sprintf('\\Lambda vs Pfa Curve for %d RMS Averages',Nrms));
grid on;


