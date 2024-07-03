% 2 nonstationary signals, 1 low pass filtered
% output together and as one reference

fs=1024;

nfft=1024;
nbuff=200;
npts=nbuff*nfft;

f=fs/nfft:fs/nfft:fs;

f0=300;
sigf0=20;

f1=105;
sigf1=10;

nzlev=.1;	% 20 dB snr

Nbar=10000;
rand_raw0=2.*sigf0.*sqrt(Nbar).*randn([1 npts]);
rand_smooth0=zeros(size(rand_raw0));
rand_raw1=2.*sigf1.*sqrt(Nbar).*randn([1 npts]);
rand_smooth1=zeros(size(rand_raw1));

alpha=(Nbar-1)/Nbar;
beta = (1 - alpha);

for n=2:npts
   rand_smooth0(n)=alpha*rand_smooth0(n-1) + beta*rand_raw0(n);
   rand_smooth1(n)=alpha*rand_smooth1(n-1) + beta*rand_raw1(n);
end;

ftrk0=zeros([1 npts]);
ftrk0=f0 + rand_smooth0;
ftrk1=zeros([1 npts]);
ftrk1=f1 + rand_smooth1;

ftrkint0=zeros(size(ftrk0));
sum0=0;
ftrkint1=zeros(size(ftrk1));
sum1=0;
for n=1:npts,
   sum0=sum0+ftrk0(n)-f0;
   ftrkint0(n)=sum0;
   sum1=sum1+ftrk1(n)-f1;
   ftrkint1(n)=sum1;
end;


nt=1:npts;
omega0=2.*pi.*(f0./fs + ftrkint0./(nt.*fs));
omega1=2.*pi.*(f1./fs + ftrkint1./(nt.*fs));

y0=sin(omega0.*nt);
y1=sin(omega1.*nt);

y1 = y1 + nzlev.*randn(size(y0));	% add noise to reference

% low pass filter y1 at 128Hz (fs/8)
m=8;
alpha=(m-1)/m;
beta=1/m;
y1lp(1)=y1(1);

%for n=2:npts,
 %  y1lp(n)=alpha*y1lp(n-1)+beta*y1(n);
%end;

y2f = y1 + y0;

[B ff tt]=specgram(y2f,nfft,fs,[],0);

figure(1);
plot(nt./fs,ftrk0,nt./fs,ftrk1);
axis([0 npts/fs 0 fs/2]);
view(90,90);


figure(2);
imagesc(ff,tt,20.*log10(abs(B')));

save y2f;
save y1;