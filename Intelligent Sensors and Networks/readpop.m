% read pop.wav

[y, fs, bits]=wavread('pop3.wav');

sound(y,fs);

[npts temp]=size(y);

n=1:npts;
t=n./fs;

step=zeros(size(y));
ramp=zeros(size(y));
quad=zeros(size(y));
step(1)=y(1);
ramp(1)=step(1);
quad(1)=ramp(1);
for nn=2:npts,
   step(nn)=step(nn-1)+y(nn);
   ramp(nn)=ramp(nn-1)+step(nn);
   quad(nn)=quad(nn-1)+ramp(nn);
end;

figure(1);
nstart=500;
nstop=1500;
subplot(4,1,1),plot(t(nstart:nstop),y(nstart:nstop),'k');
title('Pop Sound');
ylabel('Sound');
subplot(4,1,2),plot(t(nstart:nstop),step(nstart:nstop),'k');
ylabel('1st Integ');
subplot(4,1,3),plot(t(nstart:nstop),ramp(nstart:nstop),'k');
ylabel('2nd Integ');
subplot(4,1,4),plot(t(nstart:nstop),quad(nstart:nstop),'k');
xlabel('Sec');
ylabel('3rd Integ');
