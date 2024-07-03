% wavelet filters for 1/3 octave bands

fs=44100;
N=16384;    % FFT size
%N=4096;
%N=1024;

Q=3;        % Q=3 for 1/3 octave Q=1 for octave
fc=4096;
f2=fc*2^(1/(2*Q));
f1=fc/2^(1/(2*Q));

% wavelet parameters
df=.5*(f2-f1);
fa=.5*(f2+f1);
n=1.01:N;       % adding .01 avoid a log of zero
ampscale=1/fs;  % this makes the passband 0 dB with a regular fft
j=sqrt(-1);
h=ampscale.*df.*(1-cos(2.*pi.*n./N)).*exp(j.*2.*pi.*(fa/fs).*(n-N/2)).*sin(2.*pi.*(df/fs).*(n-N/2))./(2.*pi.*(df/fs).*(n-N/2));

H=fft(h,N);     % for ploting filter frequency response
F=fs.*n./N;

% lets do a check - note phase not very accurate due to df
Amp=4.2
phz=0.5
y=Amp.*(cos(2.*pi.*(fc/fs).*n + phz)+j.*sin(2.*pi.*(fc/fs).*n + phz));
outbin=h*y';
outmag=abs(outbin)
outphz=angle(outbin)

figure(1);
subplot(2,1,1),plot(n,real(h));
ylabel('h[n]');
xlabel('Sample n');
subplot(2,1,2),semilogx(F(1:N/2),20.*log10(abs(H(1:N/2))));
v=axis;
axis([20 20000 -110 10]);
ylabel('dB re 1.00');
xlabel('Hz');

figure(2);
subplot(2,1,1),plot(n,real(h));
ylabel('h[n]');
xlabel('Sample n');
subplot(2,1,2),plot(F(1:N/2),20.*log10(abs(H(1:N/2))));
v=axis;
axis([20 20000 -110 10]);
ylabel('dB re 1.00');
xlabel('Hz');

