% sampled impulse response comparison

sig=-10;
f0=25;
wd=sqrt((2*pi*f0)^2 - sig^2);

%sample rates (1000 is for true response)
fs1000=1000;
fs300=300;
fs125=125;
fs57=57;

n1000=0:500;        % these half-size arrays are for the plots
n300=0:150;
n125=0:62;
n57=0:28;

ttrue=n1000./fs1000;    % this is now a vector of time samples
ytrue=exp(sig.*ttrue).*sin(wd.*ttrue)./wd;

t57=n57./fs57;
y57=exp(sig.*t57).*sin(wd.*t57)./(exp(sig./fs57)*sin(wd/fs57));

t125=n125./fs125;
y125=exp(sig.*t125).*sin(wd.*t125)./(exp(sig./fs125)*sin(wd/fs125));

t300=n300./fs300;
y300=exp(sig.*t300).*sin(wd.*t300)./(exp(sig./fs300)*sin(wd/fs300));

figure(1);

% note, we multiply the true impulse response by wd to get them close to
% the same scale for comparison plotting
subplot(3,1,1),plot(ttrue,ytrue.*wd,'k',t57,y57,'k.');
axis([0 .5 -2 2]);
legend('true response * \omega_d','57 samples/sec');
title('Unscaled Impulse Responses');
ylabel('Response');

subplot(3,1,2),plot(ttrue,ytrue.*wd,'k',t125,y125,'k.');
axis([0 .5 -2 2]);
ylabel('Response');
legend('true response * \omega_d','125 samples/sec');

subplot(3,1,3),plot(ttrue,ytrue.*wd,'k',t300,y300,'k.');
axis([0 .5 -2 2]);
ylabel('Response');
legend('true response * \omega_d','300 samples/sec');
xlabel('Seconds');

% now for the scaled responses eqn (2.3.7)

ys57=( exp(sig/fs57)*sin(wd/fs57)/wd ).*y57;
ys125=( exp(sig/fs125)*sin(wd/fs125)/wd ).*y125;
ys300=( exp(sig/fs300)*sin(wd/fs300)/wd ).*y300;

figure(2);

subplot(3,1,1),plot(ttrue,ytrue,'k',t57,ys57,'k.');
axis([0 .5 -.01 .01]);
legend('true response','57 samples/sec');
title('Scaled Impulse Responses');
ylabel('Response');

subplot(3,1,2),plot(ttrue,ytrue,'k',t125,ys125,'k.');
axis([0 .5 -.01 .01]);
ylabel('Response');
legend('true response','125 samples/sec');

subplot(3,1,3),plot(ttrue,ytrue,'k',t300,ys300,'k.');
axis([0 .5 -.01 .01]);
ylabel('Response');
legend('true response','300 samples/sec');
xlabel('Seconds');
