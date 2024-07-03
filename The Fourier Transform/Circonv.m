% circular correlation correction via zero padding

N=256;
XN=zeros(1,N);  
YN=zeros(1,N);
YNZP=zeros(1,N);

ZN=zeros(1,N);  % circular correlation from spectral product
ZNZP=zeros(1,N);    % corrected cir cor from zero padded spectral product

RN=zeros(1,N);  % linear correlation

% Fourier transforms
FXN=zeros(1,N);
FYN=zeros(1,N);
FYNZP=zeros(1,N);
FZN=zeros(1,N);
FZNZP=zeros(1,N);

fs=1024;        % sample rate      
T=1/fs;
f=46.;
w=2*pi*f;
h0=.75;
delay=1;

% our "system is a simple delay with amplitude h0
for n=1:N,	% oldest to newest in old-school loop
	XN(n)=cos(w*n*T);             
	YN(n)=h0*cos(w*T*(n-delay)); 
	YNZP(n)=YN(n);
    % RN is the correct cross correlation
    RN(n)=h0.*cos(w.*T.*(n+delay-1))./2; % n-1 here to include 0!
end
for n=N./2:N,
	YNZP(n)=0.; % zero n/2 most recent (zero padding)
end
FXN=fft(XN);                % compute FFT's for spectral product
FYN=fft(YN);
FYNZP=fft(YNZP);
FZN=FXN.*conj(FYN);         % FZ = FX * Conj(FY) 
FZNZP=FXN.*conj(FYNZP);
ZN=ifft(FZN)./(.5.*N);      % ZN is the circular correlation
ZNZP=ifft(FZNZP)./(.5.*N);  % ZNZP is the corrected circular correlation

figure(1);
plot(RN,'k');
hold on;
plot(real(ZN),'k:o');
hold off;
axis([0 256 -.8 .8]);
xlabel(' Sample n');
ylabel('Cross Correlation');
legend('Linear XCor','Circular XCor');

figure(2);
plot(RN,'k');
hold on;
plot(real(ZNZP),'k:o');
hold off;
axis([0 256 -.8 .8]);
xlabel(' Sample n');
ylabel('Cross Correlation');
legend('Linear XCor','Corrected Cir XCor');
