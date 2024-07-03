% line array beam-steering example

thetad = 1:360;             % bearing in degrees
theta = pi.*thetad./180;    % radians

M=16;
doverl=.5;

Aw90=hanning(16);           % steering vectors
Aw60=hanning(16);
Aw30=hanning(16);

scale = 0;
for m=1:16,
 scale = scale + Aw90(m);
end;
Aw90 = Aw90.*16./scale;			% amplitude scaled
Aw60 = Aw60.*16./scale;
Aw30 = Aw30.*16./scale;

% these phase increments are added to each element phase for steering
steer60 = -j*2*pi*doverl*cos(pi*60/180);
steer30 = -j*2*pi*doverl*cos(pi*30/180);
for m=1:16,
  Aw60(m)=Aw60(m).*exp(steer60.*(m-1));
  Aw30(m)=Aw30(m).*exp(steer30.*(m-1));
end;					% steering vectors

resp90 = zeros(size(theta));
resp60 = zeros(size(theta));
resp30 = zeros(size(theta));

% simple DFT sum for beam patterns
for n=1:360,
  for m=1:16
     resp90(n) = resp90(n) + Aw90(m).*exp(j.*2.*pi.*(m-1).*doverl.*cos(theta(n)));
     resp60(n) = resp60(n) + Aw60(m).*exp(j.*2.*pi.*(m-1).*doverl.*cos(theta(n)));
     resp30(n) = resp30(n) + Aw30(m).*exp(j.*2.*pi.*(m-1).*doverl.*cos(theta(n)));
  end;
end;

figure(1);
plot(abs(resp90(1:180)),'k');
hold on;
plot(abs(resp60(1:180)),'k--');
plot(abs(resp30(1:180)),'k-.');
hold off;

axis([0,180,0,16]);
legend('90 degrees','60 degrees','30 degrees');
ylabel('|s|');
xlabel('theta');
title('16-Element Windowed Line Array Output');
text(120,14,'d/\lambda = 0.5');

figure(2);

polar(theta,abs(resp90),'k');
hold on;
polar(theta,abs(resp60),'k--');
polar(theta,abs(resp30),'k-.');
hold off;
title('16-Element Windowed Line Array Output');
legend('90 degrees','60 degrees','30 degrees');
