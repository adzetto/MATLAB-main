% check of Bayes rule
%	feature x, class wk,  cap P is probability, lc p is density

%	P(wk|x) = p(x|wk)P(wk)
%				______________
%					p(x)

xmax=20;
npts=100;

sig1=1.25;		% training std for class 1
m1=-2;			% training mean for class 1

sig2=1.;		% training std for class 2
m2=0;			% training mean for class 2

sig3=1.5;		% training std for class 3
m3=3;			% training mean for class 3

% entire class set
msum=(m1+m2+m3)/3;
sigsum=sqrt(sig1^2+sig2^2+sig3^2);


x1=-2;
x2=1;
x3=15;
xs=x3;			% our data sample

pw1=1/3;
pw2=1/3;
pw3=1/3;

x=-xmax:xmax/(npts-1):xmax;

pxw1=(1/(sqrt(2*pi)*sig1)).*exp(-.5.*((x-m1)./sig1^2).^2);
pxw2=(1/(sqrt(2*pi)*sig2)).*exp(-.5.*((x-m2)./sig2^2).^2);
pxw3=(1/(sqrt(2*pi)*sig3)).*exp(-.5.*((x-m3)./sig3^2).^2);

pxsum=(1/(sqrt(2*pi)*sigsum)).*exp(-.5.*((x-msum)./sigsum^2).^2);

figure(1);
plot(x,pxw1,'k:',x,pxw2,'k--',x,pxw3,'k-.',x,pxsum,'k');
hold on;
plot(x1,0,'k*',x2,0,'k^',x3,0,'ks');
hold off;
legend('p[x|w_1]','p[x|w_2]','p[x|w_3]','All Training','Data A','Data B','Data C');
ylabel('Density');
xlabel('Feature Value');

pxw1s=(1/(sqrt(2*pi)*sig1))*exp(-.5*((xs-m1)/sig1^2)^2)*pw1;
pxw2s=(1/(sqrt(2*pi)*sig2))*exp(-.5*((xs-m2)/sig2^2)^2)*pw2;
pxw3s=(1/(sqrt(2*pi)*sig3))*exp(-.5*((xs-m3)/sig3^2)^2)*pw3;
pxsums=(1/(sqrt(2*pi)*sigsum))*exp(-.5*((xs-msum)/sigsum^2).^2);

pxwsum=pxw1s+pxw2s+pxw3s;

pclass1=pxw1s/pxwsum;
pclass2=pxw2s/pxwsum;
pclass3=pxw3s/pxwsum;

pclsnot=pxsums*(sqrt(2*pi)*sigsum);	% normalize to 1 at mean
nstd=sqrt(-2*log(pclsnot));

d1=((xs-m1)/sig1^2)^2;
d2=((xs-m2)/sig2^2)^2;
d3=((xs-m3)/sig3^2)^2;
dsum=((xs-msum)/sigsum^2)^2;

[pclass1 pclass2 pclass3 pclsnot;
   d1 d2 d3 dsum ]'
nstd