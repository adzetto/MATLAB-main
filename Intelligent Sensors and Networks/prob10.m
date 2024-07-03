% prob10

t=1:60;
tf=20;
sfosr = 28.7;

% survivor rate
R=.5.*ones(size(t)) - (1/pi).*atan((t-20.*ones(size(t)))./sfosr);

% hazard rate
ot=ones(size(t));

H=2.*sfosr.*((pi.*ot-2.*atan((t-20.*ot)./sfosr)).^-1).*((t-20.*ot).^2 + (sfosr^2).*ot).^-1;

figure(1);

[ax,h1,h2]=plotyy(t,R,t,H);
set(get(ax(1),'Ylabel'),'String','R(t)');
set(h1,'Color','k');
set(ax(1),'YColor','k');
set(get(ax(2),'Ylabel'),'String','H(t)');
set(h2,'Color','k');
set(ax(2),'YColor','k');
xlabel('T min');
title('Survivor R(t) and Hazard H(t) Rates');