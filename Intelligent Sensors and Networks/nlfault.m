% nonlinear fault trajectory model

t=1:200;


x1=20 + .02.*(70.*t + .2.*t.^2 - .014.*t.^3 + .000075.*t.^4);
x2=20 + .02.*(50.*t + .5.*t.^2);
x3=20 + 40.*(1. - exp(-t./10));

xdata=x1 + 10.*randn(size(x1));% std of 10

FP=fopen('prognostics.dat','wt');
for nn=1:200,
   fprintf(FP,'%e %e \n',t(nn),xdata(nn));
end;
fclose(FP);

figure(1);
plot(t,x1,'k:',t,x2,'k--',t,x3,'k');
ylabel('Deg C');
xlabel('Hours');
axis([0 200 0 400]);
legend('Regulator Failure','Resistor Failure','Healthy Circuit');

% run kptrack2

npts=40;

xkp40=72.0027;	% from kptrack2
vkp40=0.0209;	% deg/min
akp40=-1.8238e-7;

tres=round(vkp40*60/.02);
res_tail=xkp40+0.01.*([tres:tres+200-npts+1]).^2 - 0.01*tres^2;
reg_tail=xkp40+0.0000015.*([tres:tres+200-npts+1]).^4 - 0.0000015*tres^4;

pred_tail(1)=xkp40;
for nn = npts+2:200,
   pred_tail(nn-npts)=pred_tail(nn-npts-1)+vkp40*60;
   vkp40=vkp40 + akp40*60;
end;


figure(2);
plot(t(1:npts),xdata(1:npts),'k+');
axis([0 200 0 400]);
hold on;

text(25, 310, 'Failure');
text(25, 160, 'Repair');

% only if kptrack2 has been run...
plot(t(1:npts),xkp(1,1:npts),'k');

plot(t(npts+1:200),res_tail(1:200-npts),'k--');
plot(t(npts+1:200),reg_tail(1:200-npts),'k:');
plot(t(npts+1:200),pred_tail(1:200-npts),'k-.');


line([0 200], [300 300]);
line([0 200], [150 150]);

hold off;
ylabel('Deg C');
xlabel('Hours');

legend('Measurements','Tracked Data','Resistor Fault','Regulator Fault','Linear Predicted');


npts=140;

xkp140=89.5909;	% from kptrack2
vkp140=0.0147;	% deg/min
akp140=3.2236e-7;

tres=round(vkp140*60/.02);	% tag on tails w/same slope
res_tail=xkp140+0.01.*([tres:tres+200-npts+1]).^2 - 0.01*tres^2;
reg_tail=xkp140+0.0000015.*([tres:tres+200-npts+1]).^4 - 0.0000015*tres^4;

pred_tail(1)=xkp140;
for nn = npts+2:200,
   pred_tail(nn-npts)=pred_tail(nn-npts-1)+vkp140*60;
   vkp140=vkp140 + akp140*60;
end;


figure(3);
plot(t(1:npts),xdata(1:npts),'k+');
axis([0 200 0 400]);
hold on;

text(25, 310, 'Failure');
text(25, 160, 'Repair');

% only if kptrack2 has been run...
plot(t(1:npts),xkp(1,1:npts),'k');

plot(t(npts+1:200),res_tail(1:200-npts),'k--');
plot(t(npts+1:200),reg_tail(1:200-npts),'k:');
plot(t(npts+1:200),pred_tail(1:200-npts),'k-.');


line([0 200], [300 300]);
line([0 200], [150 150]);

hold off;
ylabel('Deg C');
xlabel('Hours');

legend('Measurements','Tracked Data','Resistor Fault','Regulator Fault','Linear Predicted');



npts=150;

xkp150=129.5701;	% from kptrack2
vkp150=0.0440;	% deg/min
akp150=1.6602e-7;

tres=round(vkp150*60/.01);
res_tail=xkp150+0.01.*([tres:tres+200-npts+1]).^2 - 0.01*tres^2;
tres=110;
reg_tail=xkp150+0.0000015.*([tres:tres+200-npts+1]).^4 - 0.0000015*tres^4;

pred_tail(1)=xkp150;
for nn = npts+2:200,
   pred_tail(nn-npts)=pred_tail(nn-npts-1)+vkp150*60;
   vkp150=vkp150 + akp150*60;
end;


figure(4);
plot(t(1:npts),xdata(1:npts),'k+');
axis([0 200 0 400]);
hold on;

text(25, 310, 'Failure');
text(25, 160, 'Repair');

% only if kptrack2 has been run...
plot(t(1:npts),xkp(1,1:npts),'k');

plot(t(npts+1:200),res_tail(1:200-npts),'k--');
plot(t(npts+1:200),reg_tail(1:200-npts),'k:');
plot(t(npts+1:200),pred_tail(1:200-npts),'k-.');


line([0 200], [300 300]);
line([0 200], [150 150]);

hold off;
ylabel('Deg C');
xlabel('Hours');

legend('Measurements','Tracked Data','Resistor Fault','Regulator Fault','Linear Predicted');


