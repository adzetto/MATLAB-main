m=input('mean? ');
sig=input('sigma? ');
sigsq=sig.^2;

npts=300;
nsig=5;
n=1:npts;
x=2.*nsig.*sig.*(n-npts/2)./npts;
gpdf=(1./(sig.*sqrt(2.*pi))).*exp(-.5.*((x-m).^2)./sigsq);

figure(1);
plot(x,gpdf,'k');
ylabel('p_Y(x)');
xlabel('x');
titbuf=sprintf('Gaussian PDF \\sigma_Y =%.2f m_Y=%.2f',sig,m); % note the TeX
title(titbuf);

nsig=2.*nsig;
n=1:npts;
xc=(n./npts).*nsig;
cs1=(1./sqrt(2.*pi)).*xc.^(-.5).*exp(-.5.*xc);
cs2=.5*exp(-.5.*xc);
m=4;
cs4=(m./((2.^m).*gamma(m).*(sigsq.^m))).*(m.*xc).^(m-1).*exp(-.5.*(m.*xc./sigsq));
m=8;
cs8=(m./((2.^m).*gamma(m).*(sigsq.^m))).*(m.*xc).^(m-1).*exp(-.5.*(m.*xc./sigsq));
m=64;
cs64=(m./((2.^m).*gamma(m).*(sigsq.^m))).*(m.*xc).^(m-1).*exp(-.5.*(m.*xc./sigsq));

figure(2);
plot(xc,cs1,'k');
hold on;
plot(xc,cs2,'k:');
plot(xc(1:6:npts),cs4(1:6:npts),'k-^');
plot(xc,cs8,'k-.');
plot(xc(1:3:npts),cs64(1:3:npts),'k-o');
hold off;
axis([0 5 0 2]);
ylabel('p(z|v)');
xlabel('z=x^2');
legend('v=1','v=2, M=1','v=4, M=2','v=8, M=4','v=64, M=32');

sumcs1(1)=cs1(1);
sumcs2(1)=cs2(1);
sumcs4(1)=cs4(1);
sumcs8(1)=cs8(1);
sumcs64(1)=cs64(1);
for n=2:npts,
	sumcs1(n)=sumcs1(n-1)+cs1(n);
	sumcs2(n)=sumcs2(n-1)+cs2(n);
	sumcs4(n)=sumcs4(n-1)+cs4(n);
	sumcs8(n)=sumcs8(n-1)+cs8(n);
	sumcs64(n)=sumcs64(n-1)+cs64(n);
end

sumcs1=sumcs1.*nsig/npts;
sumcs2=sumcs2.*nsig./npts;
sumcs4=sumcs4.*nsig./npts;
sumcs8=sumcs8.*nsig./npts;
sumcs64=sumcs64.*nsig./npts;

figure(3);
plot(xc,sumcs1,'k');
hold on;
plot(xc,sumcs2,'k:');
plot(xc,sumcs4,'k-.');
plot(xc(1:6:npts),sumcs8(1:6:npts),'k-o');
plot(xc(1:6:npts),sumcs64(1:6:npts),'k-*');
hold off;

ylabel('P(x^2|v)');
xlabel('z=x^2');
legend('v=1','v=2, M=1','v=4, M=2','v=8, M=4','v=64, M=32');
axis([0 5 0 1.1]);

figure(4);
xclog=20.*log10(xc)./(2.*sig^2);
plot(xclog,sumcs1,'k');
hold on;
plot(xclog,sumcs2,'k:');
plot(xclog,sumcs4,'k-.');
plot(xclog(1:6:npts),sumcs8(1:6:npts),'k-o');
plot(xclog(1:6:npts),sumcs64(1:6:npts),'k-*');
hold off;

ylabel('P(x^2|v)');
xlabel('20log_{10}(x^2/2{\sigma}^2)');
legend('v=1','v=2, M=1','v=4, M=2','v=8, M=4','v=64, M=32');
axis([-10 6 0 1.1]);


