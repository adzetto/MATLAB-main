% lets plot 100 points of a damped sine wave

npts=500;
n = 0:npts;      % 1 x npts+1 vector going from 0 to npts

pi=acos(-1);    % old school - pi built in in Matlab now
j=sqrt(-1);     % old school - j built into Matlab now

sig=-10.        % damping in Nepers
f=25;           % 10 Hz
omega=2*pi*f;   % 20 pi radians per second
fs=1000;        % our sample rate

T=n./fs;        % actual time for each sample
                % note: ./ is a element-by element divide of vector n

y = exp( (sig + j*omega).*T );  % y is complex 1 x npts+1 vector now
figure(1);
plot(T,real(y),'k');
ylabel('Response');
xlabel('Seconds');
title('25 Hz -10 Nepers Sinusoid');