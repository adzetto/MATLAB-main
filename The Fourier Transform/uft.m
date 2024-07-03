npts=input('Enter npts:');
std=input('Enter std:');
fofs=input('Enter f/fs:');
t=zeros([1,npts]);
f=zeros([1,npts]);
dt=zeros([1,npts]);
yf=zeros([1,npts]);
y=zeros([1,npts]);
lomb=zeros([1,npts]);

rnd=randn([1,npts]);

t=1:npts;
t=t+std.*rnd;

%f0=128;
%fs=1024;
fs=1;
%ang=2*pi*f0/fs;
ang=2*pi*fofs;
for n=1:npts;
   y(n)=sin(14+ang*t(n));
end

figure(1);
plot(t,y,'k',t,y,'ko');
xlabel('Sample n');
ylabel('y[n]');

j=sqrt(-1);
for m=1:npts,
   f(m)=m*fs/npts;
   ang=-j*2*pi*f(m)/fs;
   ytemp=0.;
   for n=1:npts,
       ytemp=ytemp+y(n)*exp(ang*t(n));
   end
   yf(m)=ytemp;	
end

yf=(sqrt(2./npts)).*yf;
figure(2);
plot(f,20.*log10(abs(yf)),'k'); % UFT
xlabel('f/f_s');
ylabel('UFT dB');

yfft=(1./sqrt(npts)).*fft(y);
figure(3);
plot(f,20.*log10(abs(yfft)),'k'); % FFT
xlabel('f/f_s');
ylabel('FFT dB');

% Lomb Periodogram

ybar=mean(y);
ys=0;
for n=1:npts,
  ys=ys+(y(n)-ybar)^2;
end
ys=ys./(npts-1);

for m=1:npts,
  w=m.*2.*pi./npts;
  sh=0;
  ch=0;
  s2=0;
  c2=0;
  for n=1:npts,
    ang=w*t(n);
    ang2=2*ang;
    sh=sh+(y(n)-ybar)*sin(ang);
    ch=ch+(y(n)-ybar)*cos(ang);
    s2=s2+sin(ang2);
    c2=c2+cos(ang2);
  end
  tau=(.5/w)*atan2(s2,c2);
  c3=c2*cos(2*w*tau);
  s3=s2*sin(2*w*tau);	
  a=1./sqrt(.5*(npts+c3+s3));
  b=1./sqrt(.5*(npts-c3-s3));
  wt=w.*tau;
  lomb(m)=a*(ch*cos(wt)+sh*sin(wt))+j*b*(sh*cos(wt)-ch*sin(wt));
end
lomb=(.5/ys)*lomb;

figure(4);
plot(f,20.*log10(abs(lomb)),'k');   %Lomb
hold on;
plot(f,20.*log10(abs(yf)),'k:');    %UFT
hold off;
xlabel('f/f_s');
ylabel('dB');
legend('Lomb','UFT');

