% read eeee and oooo wave files

[eeee fs]=wavread('eeee.wav',8192);
[oooo fs]=wavread('oooo.wav',8192);

eeeeF=20.*log(abs(fft(eeee)));
ooooF=20.*log(abs(fft(oooo)));

f=0:fs/8192:fs-fs/8192;

enpk=0;
epklev=zeros([200 1]);
epkfrq=zeros([200 1]);
for n=20:4000;
   if (eeeeF(n)>eeeeF(n-1)) & (eeeeF(n) > eeeeF(n+1))
      bk=(sum(eeeeF(n-8:n+8))-sum(eeeeF(n-1:n+1)))/14;
      if (eeeeF(n) > (bk + 20)) & enpk <= 200
         enpk=enpk+1;
         epklev(enpk)=eeeeF(n);
         epkfrq(enpk)=f(n);
      end;
   end;
end;

onpk=0;
opklev=zeros([200 1]);
opkfrq=zeros([200 1]);
for n=20:4000;
   if (ooooF(n)>ooooF(n-1)) & (ooooF(n) > ooooF(n+1))
      bk=(sum(ooooF(n-8:n+8))-sum(ooooF(n-1:n+1)))/14;
      if (ooooF(n) > (bk + 20)) & onpk <= 200
         onpk=onpk+1;
         opklev(onpk)=ooooF(n);
         opkfrq(onpk)=f(n);
      end;
   end;
end;
            

figure(1);
subplot(2,1,1),plot(f(1:4096),eeeeF(1:4096),'k',epkfrq(1:enpk),epklev(1:enpk),'ko');
ylabel('dB')
title('eeee Sound Spectrum');
axis([0 4000 -100 100]);
subplot(2,1,2),plot(f(1:4096),ooooF(1:4096),'k',opkfrq(1:onpk),opklev(1:onpk),'ko');
ylabel('dB');
xlabel('Hz');
title('oooo Sound Spectrum');
axis([0 4000 -100 100]);
