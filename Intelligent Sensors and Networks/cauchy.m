% Cauchy pdf comparison to Gaussian

s=10;
sfsr=10;

t=[-100:100];

% Gaussian pdf
pdfg=(1/(s*sqrt(2*pi))).*exp(-.5.*(t.^2)./(s^2));

% Cauchy pdf

pdfc=(sfsr/pi).*(t.^2+ones(size(t)).*(sfsr^2)).^-1;

figure(1);

plot(t,pdfg,'k--',t,pdfc,'k');
xlabel('z');
ylabel('pdf');
legend('Gaussian','Cauchy');