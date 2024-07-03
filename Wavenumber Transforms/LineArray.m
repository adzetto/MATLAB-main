% simple line array beam

for n=0:180,
    thetad(n) = n;
    theta=pi*thetad/180;
end;
theta = pi.*(0:180)./180;

M=16;           % number of elements
%doverl1=.05;    % various spacings d over wavelength l
%doverl2=.1;
%doverl3=.5;
doverl1=.85;
doverl2=1.;
doverl3=1.5;

A=ones([1,16]);
resp1 = zeros(size(theta));
resp2 = zeros(size(theta));
resp3 = zeros(size(theta));

% n is the range over theta, m is the array element
% resp is a simple DFT sum
for n=1:181,
  for m=1:16
     resp1(n) = resp1(n) + exp(j*2*pi*(m-1)*doverl1*cos(theta(n)));
     resp2(n) = resp2(n) + exp(j*2*pi*(m-1)*doverl2*cos(theta(n)));
     resp3(n) = resp3(n) + exp(j*2*pi*(m-1)*doverl3*cos(theta(n)));
  end;
end;

% here's another way to sum the array elements using the vector features of
% Matlab
resp1 = zeros(size(theta));
resp2 = zeros(size(theta));
resp3 = zeros(size(theta));
for m=1:16
     resp1 = resp1 + exp((j*2*pi*(m-1)*doverl1).*cos(theta));
     resp2 = resp2 + exp((j*2*pi*(m-1)*doverl2).*cos(theta));
     resp3 = resp3 + exp((j*2*pi*(m-1)*doverl3).*cos(theta));
end;


figure(1);
plot(abs(resp1),'k');
hold on;
plot(abs(resp2),'k--');
plot(abs(resp3),'k-.');
axis([0,180,0,16]);
%legend('d/\lambda=0.05','d/\lambda=0.10','d/\lambda=0.50');
%tst=sprintf('d/\\lambda=%.3f',doverl1);
%legend(tst,'d/\lambda=1.00','d/\lambda=1.50');
legend('d/\lambda=0.85','d/\lambda=1.00','d/\lambda=1.50');
ylabel('|s|');
xlabel('theta');
title('16-Element Line Array Output');
hold off;