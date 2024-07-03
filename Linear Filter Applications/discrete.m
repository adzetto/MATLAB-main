m=1;
k=314;
r=4;
sig=r/(2*m);
wd=sqrt(k/m-sig^2);

km=k/m;
rm=r/m;
fs=75;
%fs=500;
T=1/fs;
Tm=T/m;
Trm=T*rm;
c=m^2/k;
A=[1 T; (-T.*km) (-Trm+1)];
B=[0; Tm];
for n=1:500,
  X(:,n)=(A^(n-1))*B;
  yd(n)=(1/T).*X(1,n);
  y(n)=(1/wd)*exp(-sig.*n.*T).*sin(wd.*n.*T);
end;
t=[1:500].*T;

figure(1);
%plot(X(1,:));
plot(t,y,'k');
hold on;
plot(t,yd,'k:');
hold off;
xlabel('sec');
legend('Continuous','Discrete');
title('Force Unit Impulse Response');
