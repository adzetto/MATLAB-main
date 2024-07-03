% butterworth responses

F=1:100:2000;   % Hz
[nrow,npts]=size(F);

w=(2*pi).*F;
s=j.*w;

fc=1000;
wc=2*pi*fc;

[B1, A1]=butter(1,wc,'s');
[B2, A2]=butter(2,wc,'s');
[B4, A4]=butter(4,wc,'s');
[B8, A8]=butter(8,wc,'s');

A1(1)=1;
H1=zeros(size(F));
num=zeros(size(F));
den=zeros(size(F));
n=1;
for k=1:n+1,
    num=num+B1(k).*s.^(n+1-k);
    den=den+A1(k).*s.^(n+1-k);
end;
H1=num./den;

A2(1)=1;
H2=zeros(size(F));
num=zeros(size(F));
den=zeros(size(F));
n=2;
for k=1:n+1,
    num=num+B2(k).*s.^(n+1-k);
    den=den+A2(k).*s.^(n+1-k);
end;
H2=num./den;

A4(1)=1;
H4=zeros(size(F));
num=zeros(size(F));
den=zeros(size(F));
n=4;
for k=1:n+1,
    num=num+B4(k).*s.^(n+1-k);
    den=den+A4(k).*s.^(n+1-k);
end;
H4=num./den;

A8(1)=1;
H8=zeros(size(F));
num=zeros(size(F));
den=zeros(size(F));
n=8;
for k=1:n+1,
    num=num+B8(k).*s.^(n+1-k);
    den=den+A8(k).*s.^(n+1-k);
end;
H8=num./den;

figure(1);
subplot(2,1,1),plot(F,(180/pi).*unwrap(angle(H1)),'k',...
    F,(180/pi).*unwrap(angle(H2)),'k:',...
    F,(180/pi).*unwrap(angle(H4)),'k-o',...
    F,(180/pi).*unwrap(angle(H8)),'k-^');
legend('1st order','2nd order', '4th order', '8th order','Location','SouthWest');
title('Butterworth Filter Responses');
ylabel('Phase [deg]');
xlabel('Hz');
%subplot(2,1,2),plot(F,20.*log10(abs(H1)),F,20.*log10(abs(H2)),...
 %   F,20.*log10(abs(H4)),F,20.*log10(abs(H8)));
subplot(2,1,2),plot(F,abs(H1),'k',F,abs(H2),'k:',...
    F,abs(H4),'k-o',F,abs(H8),'k-^');
ylabel('Magnitude');
xlabel('Hz');

dw=diff(w);   % dw
gd1=-1000.*diff(unwrap(angle(H1)))./dw;   % group delay in ms
gd2=-1000.*diff(unwrap(angle(H2)))./dw;
gd4=-1000.*diff(unwrap(angle(H4)))./dw;
gd8=-1000.*diff(unwrap(angle(H8)))./dw;
Fd=F(1:npts-1);     % we loose 1 point from differential operation

pd1=-1000.*unwrap(angle(H1))./w;          % phase delay in ms
pd2=-1000.*unwrap(angle(H2))./w;
pd4=-1000.*unwrap(angle(H4))./w;
pd8=-1000.*unwrap(angle(H8))./w;

figure(2);
subplot(2,1,1),plot(Fd,gd1,'k',Fd,gd2,'k:',Fd,gd4,'k-o',Fd,gd8,'k-^');
ylabel('Group Delay[msec]');
xlabel('Hz');
title('butterworth Phase and Time Delay Responses');
legend('1st order','2nd order', '4th order', '8th order','Location','NorthEast');

subplot(2,1,2),plot(F,pd1,'k',F,pd2,'k:',F,pd4,'k-o',F,pd8,'k-^');
ylabel('Phase Delay[msec]');
xlabel('Hz');

