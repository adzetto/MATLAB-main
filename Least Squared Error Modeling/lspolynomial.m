% basic ls poly fits

xdata = [1 2 3 4]';
ydata = [1.1 6 8 26.5]';

I4 = zeros([4 4]);
I4(1)=1;
I4(6)=1;
I4(11)=1;
I4(16)=1;

% linear fit

x1 = [ xdata'; ones(size(xdata))']';
x1x1=x1'*x1;
x1inv=x1x1^(-1);
h1 = x1inv*x1'*ydata;
px1=x1*x1inv*x1';
e1=(I4-px1)*ydata;
mse1=.25.*e1'*e1;

% quadratic fit

x2 = [(xdata.^2)'; xdata'; ones(size(xdata))']';
x2x2=x2'*x2;
x2inv=x2x2^(-1);
h2 = x2inv*x2'*ydata;
px2=x2*x2inv*x2';
e2=(I4-px2)*ydata;
mse2=.25.*e2'*e2;

% cubic fit
x3 = [(xdata.^3)'; (xdata.^2)'; xdata'; ones(size(xdata))']';
x3x3=x3'*x3;
x3inv=x3x3^(-1);
h3 = x3inv*x3'*ydata;
px3=x3*x3inv*x3';
e3=(I4-px3)*ydata;
mse3=.25.*e3'*e3;

% plot data for curves

xp=linspace(xdata(1),xdata(4),50)';
xp1 = [xp'; ones(size(xp))']';
xp2 = [(xp.^2)'; xp'; ones(size(xp))']';
xp3 = [(xp.^3)'; (xp.^2)'; xp'; ones(size(xp))']';

y1=xp1*h1;
y2=xp2*h2;
y3=xp3*h3;

figure(1);

plot(xdata,ydata,'ko');
hold on;
plot(xp,y1,'k:');
plot(xp,y2,'k-.');
plot(xp,y3,'k');

hold off;
legend('observations','linear fit','quadratic fit','cubic fit');
xlabel('x');
ylabel('y');
