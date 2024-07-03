% transfer function and coherence example

fs=1024;    % sample rate
T=1/fs;
npts=1024;
n=1:npts;
f=0:fs/npts:fs*(1-1/npts);  % frequency array in Hz starts at 0 Hz

% system poles and zeros defined by frequency and damping
fz1=240.5;
sigz1=-5;

fp1=374.5;
sigp1=-5;

fp2=60.5;
sigp2=-5;

j=sqrt(-1);     % old school

% s-plane poles and zeros
sz1=sigz1+j.*2.*pi.*fz1;
sz1c=sigz1-j.*2.*pi.*fz1;

sp1=sigp1+j.*2.*pi.*fp1;
sp1c=sigp1-j.*2.*pi.*fp1;

sp2=sigp2+j.*2.*pi.*fp2;
sp2c=sigp2-j.*2.*pi.*fp2;


%
% digital domin
%

% linear mapping s to z planes
zz1=exp(sz1.*T);
zz1c=exp(sz1c.*T);
zp1=exp(sp1.*T);
zp1c=exp(sp1c.*T);
zp2=exp(sp2.*T);
zp2c=exp(sp2c.*T);

% partial fraction expansion coeffs
Az=((zp1-zz1)*(zp1-zz1c))/((zp1-zp1c)*(zp1-zp2)*(zp1-zp2c));
Bz=((zp1c-zz1)*(zp1c-zz1c))/((zp1c-zp1)*(zp1c-zp2)*(zp1c-zp2c));
Cz=((zp2-zz1)*(zp2-zz1c))/((zp2-zp1)*(zp2-zp1c)*(zp2-zp2c));
Dz=((zp2c-zz1)*(zp2c-zz1c))/((zp2c-zp1)*(zp2c-zp1c)*(zp2c-zp2));

% ta-da!  impulse response in z domain - size of FFT in n
hz=Az*zp1.^n+Bz*zp1c.^n+Cz*zp2.^n+Dz*zp2c.^n;

wind=ones([1,npts]);
winflag=input('Hanning Window (y/n): ','s');
if winflag == 'y' | winflag == 'Y'
  wind=sqrt(3/8).*(1-cos(2*pi*n/npts));    % sqr 3/8 broadband normalization!
  %wind=(1-cos(2*pi*n/npts));    % no broadband normalization!
end

Gxy=zeros([1,npts]);    % cross spectrum
Gxx=zeros([1,npts]);    % input auto spectrum
Gyy=zeros([1,npts]);    % output auto spectrum

nd=input('number of averages?: ');
insnr=input('Input noise level: ');         % this is input corrupting noise
outnlev=input('Output noise level: ');      % this is output corrupting noise

H=fft(hz);  % true transfer function - note no 1/N because of how we defined hz
U=randn(npts*nd,1);
V=filter(hz,1,U);

% nd is the number of data blocks
for k=1:nd,

%    U=randn(npts,1);                        % noiseless signals
%    Uf=fft(U);
%    Vf=H.*Uf.';
%    V=ifft(Vf).';
    
    X=wind.*(U((k-1)*npts+1:k*npts)+insnr.*randn(npts,1))';      % noisy signals
    FX=fft(X)./npts;

    Y=wind.*(V((k-1)*npts+1:k*npts)+outnlev.*randn(npts,1))';
    FY=fft(Y)./npts;

    Gxy=Gxy+2.*conj(FX).*FY./nd;    % factor of 2 applied for non-zero Hz bins
    Gxx=Gxx+2.*conj(FX).*FX./nd;    % this is actually wrong at 0 Hz but OK in Hxy
    Gyy=Gyy+2.*conj(FY).*FY./nd;

end

% given cross and autospectra its easy to get the transfer fcn and
% coherence
Hxy=Gxy./Gxx;
gammxysq=conj(Gxy).*Gxy./(Gxx.*Gyy);

errlim=1e-10;   % this is so we don't divide by zero when nd=1
snrsq=gammxysq./(ones(size(gammxysq))-gammxysq+errlim);    % see how easy snr-squared is?

F=f(1:npts/2);  % convenient plot range
rad2deg=180/pi;

titbuf=sprintf('Num ave: %d Win:%s Nzin:%.2f Nzout:%.2f',nd,winflag,insnr,outnlev);
figure(1);
subplot(3,1,1),plot(F,rad2deg.*angle(Hxy(1:npts/2)),'k:',F,rad2deg.*angle(H(1:npts/2)),'k');
axis([0 fs/2 -200 200]);
ylabel('< H_{xy}(f) deg');
title(titbuf);
subplot(3,1,2),plot(F,20.*log10(abs(Hxy(1:npts/2))),'k:',F,20.*log10(abs(H(1:npts/2))),'k');
axis([0 fs/2 -50 50]);
ylabel('|H_{xy}(f)| dB');
subplot(3,1,3),plot(F,gammxysq(1:npts/2),'k:');
axis([0 fs/2 0 1.1]);
ylabel('{\gamma_{xy}}^2(f)');
xlabel('Hz');

figure(2);
plot(F,2.*conj(FX(1:npts/2)).*FX(1:npts/2),'k:');
hold on;
plot(F,Gxx(1:npts/2),'k');
hold off;
ylabel('G_{xx}(f)');
xlabel('Hz');
v=axis;
title(titbuf);
axis([0 fs/2 0 v(4)]);
lbuff=sprintf('M=%d',nd);
legend('M=1',lbuff);

figure(3);
subplot(2,1,1),plot(F,gammxysq(1:npts/2),'k');
axis([0 fs/2 0 1.1]);
ylabel('{\gamma_{xy}}^2(f)');
title(titbuf);
subplot(2,1,2),plot(F,10.*log10(snrsq(1:npts/2)),'k');
ylabel('SNR dB');
xlabel('Hz');
axis([0 fs/2 -20 60]);


% lets show the error bounds on the transfer fcn mag and phase.
sigH=abs(Hxy(1:npts/2)).*sqrt(1./(2.*nd.*snrsq(1:npts/2)));
sigT=atan(sqrt(1./(2.*nd.*snrsq(1:npts/2))));

sigTp=rad2deg.*(angle(Hxy(1:npts/2)) + sigT);
sigTm=rad2deg.*(angle(Hxy(1:npts/2)) - sigT);

sigHp=abs(Hxy(1:npts/2))+sigH(1:npts/2);
sigHm=abs(Hxy(1:npts/2))-sigH(1:npts/2);

figure(4);
subplot(2,1,1),plot(F,rad2deg.*angle(Hxy(1:npts/2)),'k',F,sigTp,'k:',F,sigTm,'k:');
axis([0 fs/2 -200 200]);
ylabel('< H_{xy}(f) deg');
title(titbuf);
subplot(2,1,2),plot(F,20.*log10(abs(Hxy(1:npts/2))),'k',F,20.*log10(sigHp),'k:',F,20.*log10(sigHm),'k:');
axis([0 fs/2 -50 50]);
ylabel('|H_{xy}(f)| dB');
xlabel('Hz');
