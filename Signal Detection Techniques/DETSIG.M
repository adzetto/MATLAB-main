% detection signals

n=256;
x=randn([1 n]);
y=zeros([1 n]);
xa=zeros([1 n]);
z=zeros([1 n]);

fs=10000000;    % 10 MHz
f0=3500000;
delay=100;      % 10 usec
length=10;
amp=2;
nave=10;

for i=1:length,
	y(i)=amp*cos(2*pi*f0*(i)/fs);
end;

ypow=(y*y.');

for k=1:nave,

	x=randn([1 n]);
	for i=delay:delay+length,
		x(i)=x(i)+amp*cos(2*pi*f0*i/fs);
	end;
	xa=xa+x;

%z=xcorr(y,x,'biased');

	for i=1:n/2,
		for j=1:n/2,
			z(i)=z(i)+y(j)*x(j+i);
		end;
	end;
end;
xa=xa./nave;
z=z./(ypow*nave);

t=1e6.*[1:n]./fs;
figure(1);
subplot(3,1,1), plot(t,y,'k');
ylabel('Xmit');
title(sprintf('Transmit, Receive, & Cross-Correlation with %d Ave',nave));
axis([0 25 -2*amp 2*amp]);
subplot(3,1,2), plot(t,xa,'k');
ylabel('Rec');
axis([0 25 -2*amp 2*amp]);
subplot(3,1,3), plot(t,z,'k');
ylabel('Xcor');
xlabel('usec');
axis([0 25 -2 2]);

figure(2);
subplot(3,1,1), plot(t,sqrt(y.*y),'k');
ylabel('|Xmit|');
title(sprintf('Transmit, Receive & Cross-Correlation Magnitudes with %d Ave',nave));
axis([0 25 0 amp]);
subplot(3,1,2), plot(t,sqrt(xa.*xa),'k');
ylabel('|Rec|');
axis([0 25 0 2*amp]);
subplot(3,1,3), plot(t,sqrt(z.*z),'k');
ylabel('|Xcor|');
xlabel('usec');
axis([0 25 0 1.5]);



