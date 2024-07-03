% FMCW Demo

c=1500;     % wave speed

f1=1200;   % starting freq
f2=2400;   % ending freq
Tc=2;       % chirp time in sec
beta=(f2-f1)/Tc;    % chirp rate in Hz/s

fs=10000;  % sample rate
Ts=1/fs;    % sample interval in sec

npts=round(Tc*fs);  % transmit and receive buffer size
T=[1:npts].*Ts;

fc=zeros(1,npts);
for n=1:npts,
    fc(n)=f1+n*Ts*beta;
end;

% check chirp
figure(1);
plot(fc,T);
axis([0 fs/2 0 T(npts)]);


% generate transmit signal
for n=1:npts,
    x(n)=cos(2*pi*(f1+n*Ts*beta*.5)*T(n));  % note the 1/2
end;                              % integral of freq is phase  


nobj=4;
objrng=[30 150 400 750];  % distances to onjects
%nobj=5;
%objrng=[30 150 400 402.5 750];  % distances to onjects
attnobj=.5./objrng;         % 2/R spreading there and back

objdelay=2.*objrng./c;      % delay in sec there and back
objoff=round(objdelay./Ts);

% generate received signal
y=zeros(size(x));
for k=1:nobj,
    y(objoff(k):npts)=y(objoff(k):npts)+attnobj(k).*x(1:npts-objoff(k)+1);
end;
y=y+.0005.*randn(size(y));   % some noise to make it real looking

% check the transmit and received signals
figure(2);
[xS,xf,xt,xp]=spectrogram(x,256,128,256,fs); 
subplot(2,1,1),surf(xt,xf,20.*log10(abs(xS)),'EdgeColor','none');
view(90,-90);
xlabel('Sec');
ylabel('Hz');
title('Transmitted FMCW Signal');
colormap('bone');
v=caxis;
caxis([.7*v(1) v(2)]);  % darken the background some
shading interp;

%subplot(2,1,2),spectrogram(y,256,128,1024,fs); % spotty plot!
[yS,yf,yt,yp]=spectrogram(y,256,128,256,fs); 
subplot(2,1,2),surf(yt,yf,20.*log10(abs(yS)),'EdgeColor','none');
view(90,-90);
xlabel('Sec');
ylabel('Hz');
title('Received FMCW Signal');
colormap('bone');
v=caxis;
caxis([.7*v(1) v(2)]);
shading interp;

% OK let mix them
z=x.*y;

% check mixed signal
figure(3);
%spectrogram(z,256,128,1024,fs); 
[zS,zf,zt,zp]=spectrogram(z,256,128,256,fs); 
surf(zt,zf,20.*log10(abs(zS)),'EdgeColor','none');
view(90,-90);
xlabel('Sec');
ylabel('Hz');
title('Mixed FMCW Signal');
colormap('bone');
v=caxis;
caxis([.7*v(1) v(2)]);
shading interp;

% lets look in detail
han=.5.*(1-cos((2*pi*[1:16384]./16384)));
zh=z(1:16384).*han; % Hanning window reduces leakage
Rec=fft(zh,16384);

F=[1:16384].*fs./16384; % x-axis in Hz
Rng=F.*c./(2*beta);         % x-axis in meters! (dist to obj)

figure(4);
pltpts=1200;
plot(Rng(1:pltpts),20.*log10(abs(Rec(1:pltpts))),'k');
ylabel('dB');
xlabel('Range (m)');
axis([0 Rng(pltpts) -35 35]);
grid on;
%title('FMCW Range Estimation');
% let's add a second x-axis on top
ax1=gca;
ax2=axes('Position',get(ax1,'Position'),...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none','XColor','k','YColor','k');
h12=line(F(1:pltpts),20.*log10(abs(Rec(1:pltpts))),...
    'Color','k','Parent',ax2);
xlabel('Hz');
axis([0 F(pltpts) -35 35]);


% lets zoom in where we want
figure(5);
clf;
pltbeg=400;
pltend=600;
plot(Rng(pltbeg:pltend),20.*log10(abs(Rec(pltbeg:pltend))),'k');
ylabel('dB');
xlabel('Range (m)');
axis([Rng(pltbeg) Rng(pltend) -35 10]);
grid on;
%title('FMCW Range Estimation');
% let's add a second x-axis on top
ax1=gca;
ax2=axes('Position',get(ax1,'Position'),...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none','XColor','k','YColor','k');
h12=line(F(pltbeg:pltend),20.*log10(abs(Rec(pltbeg:pltend))),'Color','k','Parent',ax2);
xlabel('Hz');
axis([F(pltbeg) F(pltend) -35 10]);
