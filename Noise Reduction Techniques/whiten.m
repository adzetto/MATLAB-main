% whitening filter and ALE for y.mat

load y;

nwts=24;

[temp npts]=size(y);

err=zeros(size(y));
yh=zeros(size(y));
yale=zeros(size(y));
a=zeros([1 nwts]);

mu=.005/nwts;

for n=nwts+1:npts,
   yh(n)=-y(n-nwts:n-1)*a.';
   err(n)=y(n)-yh(n);
   emu=mu*err(n);
   a=a-emu.*y(n-nwts:n-1);
%   yale(n)=-yh(n-nwts:n-1)*a.';
end;

[Berr ff tt]=specgram(err,nfft,fs,[],0);
[Byh ff tt]=specgram(yh,nfft,fs,[],0);
[By ff tt]=specgram(y,nfft,fs,[],0);



figure(1);

subplot(1,3,1),imagesc(ff,tt,20.*log10(abs(By')));
xlabel('Hz');
ylabel('Sec');
title('Original');
caxis([-80 50]);

subplot(1,3,2), imagesc(ff,tt,20.*log10(abs(Byh')));
xlabel('Hz');
ylabel('Sec');
title('Enhanced');
caxis([-80 50]);

subplot(1,3,3),imagesc(ff,tt,20.*log10(abs(Berr')));
xlabel('Hz');
ylabel('Sec');
title('Whitened');
caxis([-80 50]);
colorbar('vert');
text(700,55,'dB');
colormap(bone);


rec1=60;
rec2=70;

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
titlebuf=sprintf('Whitening Filter Response Average T=%d:%d',rec1,rec2);
title(titlebuf);
legend('Original','Whitened');
axis([0 nfft/2 -50 50]);

figure(3);
plot(f(1:nfft/2),20.*log10(Yf),'k');
hold on;
plot(f(1:nfft/2),20.*log10(Yfe),'k:');
hold off;
ylabel('dB');
xlabel('Hz');
titlebuf=sprintf('Adaptive Line Enhancement Response Average T=%d:%d',rec1,rec2);
title(titlebuf);
legend('Original','Enhanced');
axis([0 nfft/2 -50 50]);
