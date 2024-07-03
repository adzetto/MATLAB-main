w0=20;
sig=input('sigma?: ');
wd=sqrt(w0^2-sig^2);

for n=1:400,
 wn(n)=18+4*n/400;
 w=18+4*n/400;
 z(n)=(sig^2+j*2*sig*w-w^2+wd^2)/(j*w);	% mech z
 zw(n)=(sig^2+j*2*sig*w-w^2+wd^2);      % 1/displacement
 za(n)=(sig^2+j*2*sig*w-w^2+wd^2)/(-w^2);  % 1/accel

end;

figure(1);
plot(wn,20*log10(abs(1./z)),'r');
hold on;
plot(wn,20*log10(abs(1./zw)),'g');
plot(wn,20*log10(abs(1./za)),'y');
hold off;
