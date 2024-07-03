% line array beam patterns with and without shading

thetad = 0:181;             % theta range in degrees
theta = pi.*thetad./180;    % theta in radians

M=16;                       % 16 elements
doverl=.5;                  % d/l is max before grating lobes occur

A=ones([16,1]);             % this is a 16-sample rectangular window
Aw=hanning(16);             % this is a 16 sample Hanning window
scale = 0;
for m=1:16,
 scale = scale + Aw(m);     % here we just sum up the window for scaling
end;
Aw = Aw.*16./scale;         % this re-scales to match scale of rect window
resp = zeros(size(theta));
respw = zeros(size(theta));

% simple DFT sum at wach angle produces the beam patterns
for n=1:182,
  for m=1:16
     resp(n) = resp(n) + A(m).*exp(j.*2.*pi.*(m-1).*doverl.*cos(theta(n)));
     respw(n) = respw(n) + Aw(m).*exp(j.*2.*pi.*(m-1).*doverl.*cos(theta(n)));
  end;
end;

figure(1);
plot(abs(resp),'k');
hold on;
plot(abs(respw),'k--');
hold off;

axis([0,180,0,16]);
legend('Rectangular','Hanning');
ylabel('|s|');
xlabel('theta');
title('16-Element Windowed Line Array Output');
text(120,10,'d/\lambda = 0.5');

