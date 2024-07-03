% Eigenvalue case demo

d=.05;
numels=40;
f=150;
c=345;
k=2*pi*f/c;
j=sqrt(-1);

D=zeros([numels 1]);
for n=1:numels,
	D(n)=(n-1)*d;
end;

theta1=90;
y1=exp(j.*k.*D.*cos((pi/180).*theta1));
theta2=60;
y2=exp(j.*k.*D.*cos((pi/180).*theta2));
theta3=0;
y3=exp(j.*k.*D.*cos((pi/180).*theta3));

x0=zeros(size(y1));
y0=zeros(size(y1));
y1x=real(y1);
y1y=imag(y1);
y2x=real(y2);
y2y=imag(y2);
y3x=real(y3);
y3y=imag(y3);

figure(1);
subplot(3,1,1), plot3(D,y1x,y1y,'k');
axis([D(1) D(numels) -1 1 -1 1]);
hold on;
for n=1:numels,
	h=line([D(n) D(n)],[x0(n) y1x(n)],[y0(n) y1y(n)]);
	%set(h,'LineStyle','--');
	set(h,'Color','k');
end;
h=line([D(1) D(numels)],[0 0],[0 0]);
set(h,'Color','k');
hold off;
view(10,30);
ylabel('Real[y]');
zlabel('Imag[y]');
title('150 Hz at 90 degrees bearing');
grid on;

subplot(3,1,2), plot3(D,y2x,y2y,'k');
axis([D(1) D(numels) -1 1 -1 1]);
hold on;
for n=1:numels,
	h=line([D(n) D(n)],[x0(n) y2x(n)],[y0(n) y2y(n)]);
	%set(h,'LineStyle','--');
	set(h,'Color','k');
end;
h=line([D(1) D(numels)],[0 0],[0 0]);
set(h,'Color','k');
hold off;
view(10,30);
ylabel('Real[y]');
zlabel('Imag[y]');
title('60 degrees bearing');
grid on;

subplot(3,1,3), plot3(D,y3x,y3y,'k');
axis([D(1) D(numels) -1 1 -1 1]);
hold on;
for n=1:numels,
	h=line([D(n) D(n)],[x0(n) y3x(n)],[y0(n) y3y(n)]);
%	set(h,'LineStyle','--');
	set(h,'Color','k');
end;
h=line([D(1) D(numels)],[0 0],[0 0]);
set(h,'Color','k');
hold off;
view(10,30);
xlabel('Line Array Position D');
ylabel('Real[y]');
zlabel('Imag[y]');
title('0 degrees bearing');
grid on;

% combined wave
yt=(y1+y2+y3);
ytx=real(yt);
yty=imag(yt);

figure(2)
subplot(3,1,1), plot3(D,ytx,yty,'k');
axis([D(1) D(numels) -3 3 -3 3]);
hold on;
for n=1:numels,
	h=line([D(n) D(n)],[x0(n) ytx(n)],[y0(n) yty(n)]);
%	set(h,'LineStyle','--');
	set(h,'Color','k');
end;
h=line([D(1) D(numels)],[0 0],[0 0]);
set(h,'Color','k');
hold off;
view(10,30);
ylabel('Real[y]');
zlabel('Imag[y]');
title('150 Hz Complex Multipath Response');
h=gca;
set(h,'ZTick',[-3 0 3],'YTick',[-3 0 3],'XTick',[0 .5 1 1.5]);
grid on;

subplot(3,1,2), plot(D,abs(yt),'k');
ylabel('Abs[y]');
title('Magnitude of Multipath Response');
h=gca;
set(h,'YTick',[-180 -90 0 90 180],'XTick',[0 .5 1 1.5]);

subplot(3,1,3), plot(D,(180/pi).*angle(yt),'k');
ylabel('Degrees');
xlabel('Line Array Position D');
title('Phase of Multipath Response');
axis([0 D(numels) -180 180]);
h=gca;
set(h,'YTick',[-180 -90 0 90 180],'XTick',[0 .5 1 1.5]);

