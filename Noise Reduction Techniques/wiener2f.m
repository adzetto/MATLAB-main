% whitening filter and ALE for y.mat

load y1;
load y2f;

nwts=8;
nfft=1024;

[temp npts]=size(y1);

err=zeros(size(y1));
yh=zeros(size(y1));

a=zeros([1 nwts]);

murel=.05;
pow=std(y1(1:100))^2;
mu=murel/(nwts*pow);

for n=nwts:npts,
   yh(n)=y1(n-nwts+1:n)*a.';
   err(n)=y2f(n)-yh(n);
   emu=mu*err(n);
   a=a+emu.*y1(n-nwts+1:n);
end;

[Berr ff tt]=specgram(err,nfft,fs,[],0);
[Byh ff tt]=specgram(yh,nfft,fs,[],0);
[By ff tt]=specgram(y2f,nfft,fs,[],0);



figure(1);

subplot(1,3,1),imagesc(ff,tt,20.*log10(abs(By')));
xlabel('Hz');
ylabel('Sec');
title('Original');
caxis([-30 50]);

subplot(1,3,2), imagesc(ff,tt,20.*log10(abs(Byh')));
xlabel('Hz');
ylabel('Sec');
title('Predicted Interference');
caxis([-30 50]);

subplot(1,3,3),imagesc(ff,tt,20.*log10(abs(Berr')));
xlabel('Hz');
ylabel('Sec');
title('Canceled Interference');
caxis([-30 50]);
colorbar('vert');
nbuff=npts/nfft;
text(700,nbuff/2,'dB');
colormap(bone);
set(gcf,'color','w');


rec1=40;
rec2=45;

Yf=zeros([1 512]);
Yfe=Yf;
Yfw=Yf;

for k=rec1:rec2,
   Yf=Yf+abs(B(1:nfft/2,k)');
   Yfe=Yfe+abs(Byh(1:nfft/2,k)');
   Yfw=Yfw+abs(Berr(1:nfft/2,k)');
end;
Yf=Yf./(rec2-rec1+1);
Yfe=Yfe./(rec2-rec1+1);
Yfw=Yfw./(rec2-rec1+1);

figure(2);
plot(f(1:nfft/2),20.*log10(Yf),'k');
hold on;
plot(f(1:nfft/2),20.*log10(Yfw),'k:');
hold off;
ylabel('dB');
xlabel('Hz');
titlebuf=sprintf('Cancellation Result Average T=%d:%d',rec1,rec2);
title(titlebuf);
legend('Original','Canceled');
axis([0 nfft/2 -50 50]);
