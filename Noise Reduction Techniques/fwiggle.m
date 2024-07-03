% nonstationary signal and stationary noise canceller
% state model approach

fs=1024;

nfft=1024;
nbuff=200;
npts=nbuff*nfft;

f=fs/nfft:fs/nfft:fs;

f0=350;
sigf=10;

nzlev=.01;	% 40 dB snr

Nbar=10000;
rand_raw=2.*sigf.*sqrt(Nbar).*randn([1 npts]);
rand_smooth=zeros(size(rand_raw));

alpha=(Nbar-1)/Nbar;
beta = (1 - alpha);

for n=2:npts
   rand_smooth(n)=alpha*rand_smooth(n-1) + beta*rand_raw(n);
end;

ftrk=zeros([1 npts]);
ftrk=f0 + rand_smooth;

ftrkint=zeros(size(ftrk));
sum=0;
for n=1:npts,
   sum=sum+ftrk(n)-f0;
   ftrkint(n)=sum;
end;


nt=1:npts;
omega=2.*pi.*(f0./fs + ftrkint./(nt.*fs));

f1=120.5;
omega1=2*pi*f1/fs;
y=sin(omega.*nt) + sin(omega1.*nt);

y = y + nzlev.*randn(size(y));

[B ff tt]=specgram(y,nfft,fs,[],50);

figure(1);
plot(nt./fs,ftrk);
axis([0 npts/fs 0 fs/2]);
view(90,90);


figure(2);
imagesc(ff,tt,20.*log10(abs(B')));

save y;