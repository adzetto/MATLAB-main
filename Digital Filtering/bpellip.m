% check s-to-z mapping for complex filter

fs=4000;        % sample rate for digital system
T=1/fs;

F=1:1:2000;   % Hz
[nrow,npts]=size(F);
w=(2*pi).*F;
s=j.*w;

f1=1200;
f2=1400;
f0=.5*(f2-f1);
w0=2*pi*f0;

Wp=[2*pi*f1,2*pi*f2];

Rp=0.1;     % pass band ripple
Rs=60;      % stop band attenuation

N=8;           % this results in 2N+1 coeeficients for band pass
                % works for N even

[Bs,As]=ellip(N,Rp,Rs,Wp,'bandpass','s');    % bandpass s-plane filter
[Zs,Ps,Ks]=ellip(N,Rp,Rs,Wp,'bandpass','s');   % get zeros and poles
Hs=freqs(Bs,As,w);

% specify on z-plane
Wpz=2.*Wp./(2*pi*fs);
[Bz,Az]=ellip(N,Rp,Rs,Wpz,'bandpass');    % bandpass z-plane filter
[Zz,Pz,Kz]=ellip(N,Rp,Rs,Wpz,'bandpass');   % get zeros and poles
Hz=freqz(Bz,Az,npts);


[Bsb,Asb]=bilinear(Bs,As,fs);
Hzb=freqz(Bsb,Asb,npts);        % bilinear mapping

fp=.5*(f1+f2);
[Bsbp,Asbp]=bilinear(Bs,As,fs,fp);
Hzbp=freqz(Bsbp,Asbp,npts);        % bilinear mapping w/warping

% lets just map the poles and zeros
Zsz=exp(Zs.*T);
Psz=exp(Ps.*T);

% estimate the scale factor 2*N always even, z's in conj pairs
num=1;
den=1;
for kk=1:N,
    k=2*(kk-1)+1;
    num=num*(Psz(k)-Psz(k+1))*(Zs(k)-Zs(k+1));
    den=den*(Ps(k)-Ps(k+1))*(Zsz(k)-Zsz(k+1));
end;
scale = T*num/den; 

Bsz=scale*poly(Zsz);
Asz=poly(Psz);
Hsz=freqz(Bsz,Asz,npts);

% can't put them all on one plot
figure(1);
plot(F,20.*log10(abs(Hs)),'k',F,20.*log10(abs(Hzb)),'k:',...
    F,20.*log10(abs(Hzbp)),'k--');
legend('s-plane','bi-linear','bi-lin warp to f_0');
v=axis;
axis([v(1) v(2) -80 6]);
xlabel('Hz');
ylabel('dB');

figure(2);
plot(F,20.*log10(abs(Hs)),'k',F,20.*log10(abs(Hz)),'k:',...
    F,20.*log10(abs(Hsz)),'k--');
legend('s-plane','z-plane','linear map s-to-z');
v=axis;
axis([v(1) v(2) -80 6]);
xlabel('Hz');
ylabel('dB');

figure(3);
plot(F,abs(Hs),'k',F,abs(Hzb),'k:',...
    F,abs(Hzbp),'k--');
legend('s-plane','bi-linear','bi-lin warp to f_0');
v=axis;
axis([v(1) v(2) 0 1.2]);
xlabel('Hz');
ylabel('Response');

figure(4);
plot(F,abs(Hs),'k',F,abs(Hz),'k:',...
    F,abs(Hsz),'k--');
legend('s-plane','z-plane','linear map s-to-z');
v=axis;
axis([v(1) v(2) 0 1.2]);
xlabel('Hz');
ylabel('Response');
