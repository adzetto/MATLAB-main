% pulse train generator

% Nr/(2Nw) gives number of harmonics in freq domain!
Nr=100;	% sample # pulse repeats
Nw=10;		% pulse width in samples
beta=(Nr/Nw);		% T/N is pulse width in samples
tau=.5;	% small delay
% if beta is odd -> even harmonics
% if beta is even -> odd harmonics

fs=1000;
Ts=1/fs;
npts=8192;
n=1:npts;
t=Ts*n;
f=(fs/npts)*n;

y=(1/beta).*sin(pi.*(n-tau)./Nw)./(sin(pi.*(n-tau)./Nr));

for nn=1:npts,
   if y(nn)>1
      y(nn)=1;
   end;
   if y(nn)<-1
      y(nn)=-1;
   end;
end;

      
figure(1);
subplot(2,2,1),plot(t(1:500),y(1:500),'k');
axis([0 t(500) -1 1]);
xlabel('Sec');
titbuf=sprintf('Nr=%d Nw=%d Fs=%.f',Nr,Nw,fs);
title(titbuf);

Ys=20.*log10((sqrt(2)/npts).*abs(fft(y)));

subplot(2,2,2), plot(f(1:1000),Ys(1:1000),'k');
axis([0 f(1000) -100 0]);
ylabel('dB');
xlabel('Hz');
title('Spectrum');

% do two more
Nr=100;	% sample # pulse repeats
Nw=11.1111;		% pulse width in samples
beta=(Nr/Nw);		% T/N is pulse width in samples
tau=70.1;
tau=.5;
y2=(1/beta).*sin(pi.*(n-tau)./Nw)./(sin(pi.*(n-tau)./Nr));

for nn=1:npts,
   if y(nn)>1
      y(nn)=1;
   end;
   if y(nn)<-1
      y(nn)=-1;
   end;
end;

subplot(2,2,3),plot(t(1:500),y2(1:500),'k');
axis([0 t(500) -1 1]);
xlabel('Sec');
titbuf=sprintf('Nr=%d Nw=%.2f Fs=%.f',Nr,Nw,fs);
title(titbuf);

Ys2=20.*log10((sqrt(2)/npts).*abs(fft(y2)));

subplot(2,2,4), plot(f(1:1000),Ys2(1:1000),'k');
axis([0 f(1000) -100 0]);
ylabel('dB');
xlabel('Hz');
title('Spectrum');

% some basic investigation

Ysc=fft(y2);
pkf=zeros([1 10]);
pkv=zeros([1 10]);
pkp=zeros([1 10]);
pkt=zeros([1 10]);
npk=0;
for nn=2:1000,
   if (Ys2(nn) > Ys2(nn-1)) & (Ys2(nn) > Ys2(nn+1))
      npk=npk+1;
      pkf(npk)=f(nn);
      pkv(npk)=Ys2(nn);
      pkp(npk)=unwrap(angle(Ysc(nn)));
      pkt(npk)=pkp(npk)/(2*pi*pkf(npk)); % time delay
   end;
end;

      
