% butterworth responses

F=1:50:2000;   % Hz
[nrow,npts]=size(F);

w=(2*pi).*F;
s=j.*w;

fc=1000;
wc=2*pi*fc;

% Compare 4th order filters
order=8;
[BBut, ABut]=butter(order,wc,'s');
R=1;  % 1 dB ripple in pass band
[BCheb1, ACheb1]=cheby1(order,R,wc,'s');

Rs=60; % stop band 60 dB down
e=sqrt(10^(-Rs/10)/(1+10^(-Rs/10)));
wct2=wc*(cosh((1/order)*acosh(1/e)));     % note the needed shift in cutoff
[BCheb2, ACheb2]=cheby2(order,Rs,wct2,'s');

Rs=60; % stop band 60 dB down
[BEllip, AEllip]=ellip(order,R,Rs,wc,'s');

[BBes, ABes]=besself(order,wc);

% compute frequency response on s-plane
HBut=freqs(BBut,ABut,w);
HCheb1=freqs(BCheb1,ACheb1,w);
HCheb2=freqs(BCheb2,ACheb2,w);
HEllip=freqs(BEllip,AEllip,w);
HBes=freqs(BBes,ABes,w);

figure(1);
subplot(2,1,1),plot(F,(180/pi).*unwrap(angle(HBut)),'k',...
    F,(180/pi).*unwrap(angle(HCheb1)),'k:',...
    F,(180/pi).*unwrap(angle(HCheb2)),'k-o',...
    F,(180/pi).*unwrap(angle(HEllip)),'k-^',...
    F,(180/pi).*unwrap(angle(HBes)),'k--');
legend('Butterworth','Chebyshev I', 'Chebyshev II', 'Elliptic','Bessel',...
    'Location','SouthWest');
title('8th Order Filter Responses');
ylabel('Phase [deg]');
xlabel('Hz');
subplot(2,1,2),plot(F,abs(HBut),'k',F,abs(HCheb1),'k:',...
    F,abs(HCheb2),'k-o',F,abs(HEllip),'k-^',F,abs(HBes),'k--');
v=axis;
axis([v(1) v(2) 0 1.25]);
ylabel('Magnitude');
xlabel('Hz');

dw=diff(w);   % dw
gdBut=-1000.*diff(unwrap(angle(HBut)))./dw;   % group delay in ms
gdCheb1=-1000.*diff(unwrap(angle(HCheb1)))./dw;
gdCheb2=-1000.*diff(unwrap(angle(HCheb1)))./dw;
gdEllip=-1000.*diff(unwrap(angle(HEllip)))./dw;
gdBes=-1000.*diff(unwrap(angle(HBes)))./dw;
Fd=F(1:npts-1);     % we loose 1 point from differential operation

pdBut=-1000.*unwrap(angle(HBut))./w;          % phase delay in ms
pdCheb1=-1000.*unwrap(angle(HCheb1))./w;
pdCheb2=-1000.*unwrap(angle(HCheb2))./w;
pdEllip=-1000.*unwrap(angle(HEllip))./w;
pdBes=-1000.*unwrap(angle(HBes))./w;

figure(2);
subplot(2,1,1),plot(Fd,gdBut,'k',Fd,gdCheb1,'k:',...
    Fd,gdCheb2,'k-o',Fd,gdEllip,'k-^',Fd,gdBes,'k--');
ylabel('Group Delay[msec]');
xlabel('Hz');
title('8th Order Phase and Time Delay Responses');
legend('Butterworth','Chebyshev I', 'Chebyshev II', 'Elliptic','Bessel',...
        'Location','NorthEast');
v=axis;
axis([v(1) v(2) 0 8]);    

subplot(2,1,2),plot(F,pdBut,'k',F,pdCheb1,'k:',...
    F,pdCheb2,'k-o',F,pdEllip,'k-^',F,pdBes,'k--');
ylabel('Phase Delay[msec]');
xlabel('Hz');

figure(3);

plot(F,20.*log10(abs(HBut)),'k',F,20.*log10(abs(HCheb1)),'k:',...
    F,20.*log10(abs(HCheb2)),'k-o',F,20.*log10(abs(HEllip)),'k-^',...
    F,20.*log10(abs(HBes)),'k--');
ylabel('dB');
xlabel('Hz');
title('8th Order Filter Responses');
legend('Butterworth','Chebyshev I', 'Chebyshev II', 'Elliptic','Bessel',...
        'Location','SouthWest');
v=axis;
axis([v(1) v(2) v(3) 10]);    


