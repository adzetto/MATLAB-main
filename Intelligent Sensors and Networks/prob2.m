% class A
s=[2.5 -3.2; -3.2 5.5];
m=[5; 9];

% rotation angles
theta=.5*atan2(2*s(1,2),(s(2,2)-s(1,1)));
thetaAdeg=theta*180/pi;

R=[cos(theta) sin(theta);
   -sin(theta) cos(theta)];
RA=R;
sig=R'*s*R;
sigA=sig;

nptse=50;
x1=zeros([1 nptse+1]);    % class 1 feature 1
y1=zeros([1 nptse+1]);    % class 1 feature 2

% generate ellipse
for n=1:nptse+1, 
   theta=2*pi*n/nptse;
   x1(n)=sqrt(sig(1,1))*cos(theta);
   y1(n)=sqrt(sig(2,2))*sin(theta);
end;

ellp1=R*[x1; y1];
ellp1(1,:)=ellp1(1,:)+m(1);
ellp1(2,:)=ellp1(2,:)+m(2);


figure(1);
plot(ellp1(1,:),ellp1(2,:),'k');

% class B
s=[1.75 .25; .25 6.5];
m=[2; 5];

% rotation angles
theta=.5*atan2(2*s(1,2),(s(2,2)-s(1,1)));
thetaBdeg=theta*180/pi;

R=[cos(theta) sin(theta);
   -sin(theta) cos(theta)];
RB=R;
sig=R'*s*R;
sigB=sig;

nptse=50;
x1=zeros([1 nptse+1]);    % class 1 feature 1
y1=zeros([1 nptse+1]);    % class 1 feature 2

% generate ellipse
for n=1:nptse+1, 
   theta=2*pi*n/nptse;
   x1(n)=sqrt(sig(1,1))*cos(theta);
   y1(n)=sqrt(sig(2,2))*sin(theta);
end;

ellp2=R*[x1; y1];
ellp2(1,:)=ellp2(1,:)+m(1);
ellp2(2,:)=ellp2(2,:)+m(2);


figure(2);
plot(ellp2(1,:),ellp2(2,:),'k');

figure(3)
plot(ellp1(1,:),ellp1(2,:),'k',ellp2(1,:),ellp2(2,:),'k:');
ylabel('Feature 2');
xlabel('Feature 1');
legend('Class A', 'Class B');

