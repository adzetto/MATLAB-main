% line array beam pattern looking at end-fire response

thetad = 1:360;
theta = pi.*thetad./180;

M=16;
doverl1=.05;
doverl2=.1;
doverl3=.5;

Aw1=hanning(16);
Aw2=hanning(16);
Aw3=hanning(16);

scale = 0;
for m=1:16,
 scale = scale + Aw1(m);
end;
Aw1 = Aw1.*16./scale;			% amplitude scaled
Aw2 = Aw1;
Aw3 = Aw2;

steer1 = -j*2*pi*doverl1*cos(pi*0/180);
steer2 = -j*2*pi*doverl2*cos(pi*0/180);
steer3 = -j*2*pi*doverl3*cos(pi*0/180);
for m=1:16,
  Aw1(m)=Aw1(m).*exp(steer1.*(m-1));
  Aw2(m)=Aw2(m).*exp(steer2.*(m-1));
  Aw3(m)=Aw3(m).*exp(steer3.*(m-1));
end;					% steering vectors

resp1 = zeros(size(theta));
resp2 = zeros(size(theta));
resp3 = zeros(size(theta));

for n=1:360,
  for m=1:16
     resp1(n) = resp1(n) + Aw1(m).*exp(j.*2.*pi.*(m-1).*doverl1.*cos(theta(n)));
     resp2(n) = resp2(n) + Aw2(m).*exp(j.*2.*pi.*(m-1).*doverl2.*cos(theta(n)));
     resp3(n) = resp3(n) + Aw3(m).*exp(j.*2.*pi.*(m-1).*doverl3.*cos(theta(n)));
  end;
end;

figure(1);
plot(abs(resp1(1:180)),'k');
hold on;
plot(abs(resp2(1:180)),'k--');
plot(abs(resp3(1:180)),'k-.');
hold off;

axis([0,180,0,16]);
legend('d/l=.05','d/l=0.1','d/l=0.5');
ylabel('|s|');
xlabel('theta');
title('16-Element Windowed Line Array Output');
text(80,14,'Steering Angle 0 Degrees');

figure(2);

polar(theta,abs(resp1),'k');
hold on;
polar(theta,abs(resp2),'k--');
polar(theta,abs(resp3),'k-.');
hold off;
title('16-Element Windowed Line Array Output');
legend('d/\lambda=.05','d/\lambda=0.1','d/\lambda=0.5');
