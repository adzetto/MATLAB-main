% Neural Net classifier model

npts=1500;
m11=4.5;
s11=1;
m21=6;
s21=2.5;

m12=8;
s12=.5;
m22=10;
s22=1.25;

f11=s11.*randn([1 npts]);
f21=s21.*randn([1 npts]);
data1=[f11; f21];

f12=s12.*randn([1 npts]);
f22=s22.*randn([1 npts]);
data2=[f12; f22];

tilt1=(pi/180)*(-80);
tilt2=(pi/180)*(-20);

R1=[cos(tilt1) -sin(tilt1);
   sin(tilt1) cos(tilt1)];
R2=[cos(tilt2) -sin(tilt2);
   sin(tilt2) cos(tilt2)];

datap1=R1*data1;
datap2=R2*data2;

datap1(1,:) = datap1(1,:) + m11;
datap1(2,:) = datap1(2,:) + m21;
datap2(1,:) = datap2(1,:) + m12;
datap2(2,:) = datap2(2,:) + m22;

nptse=50;
x1=zeros([1 nptse]);
y1=zeros([1 nptse]);
x2=zeros([1 nptse]);
y2=zeros([1 nptse]);

% generate ellipse
for n=1:nptse, 
   theta=2*pi*n/nptse;
   x1(n)=s11*cos(theta);
   y1(n)=s21*sin(theta);
   x2(n)=s12*cos(theta);
   y2(n)=s22*sin(theta);
end;

ellp1=R1*[x1; y1];
ellp1(1,:)=ellp1(1,:)+m11;
ellp1(2,:)=ellp1(2,:)+m21;
ellp2=R2*[x2; y2];
ellp2(1,:)=ellp2(1,:)+m12;
ellp2(2,:)=ellp2(2,:)+m22;

figure(1);
plot(datap1(1,1:10:npts),datap1(2,1:10:npts),'kx',datap2(1,1:10:npts),datap2(2,1:10:npts),'ko');
axis([0 15 0 15]);

hold on;
h=plot(ellp1(1,:),ellp1(2,:),'k',ellp2(1,:),ellp2(2,:),'k:');
hold off;
legend(h,'x - class 1','o - class 2');
xlabel('x_1');
ylabel('x_2');

% simple neural net

a11=[.01 -.02];
a12=[.2 .025];
a13=[.5 .11];

a21=[.02 .03 .04];
a22=[.023 .043 .011];

mu=.1;

for n=1:npts,
 for q=1:1,  
   %xin1=datap1(1,n)/15;
   %xin2=datap1(2,n)/15;
   xin1=m11;
   xin2=m21;
   
   net11=a11(1)*xin1+a11(2)*xin2;	% hidden layer
   o11=1/(1+exp(-net11));
   %o11=net11;
   net12=a12(1)*xin1+a12(2)*xin2;
   o12=1/(1+exp(-net12));
   %o12=net12;
   net13=a13(1)*xin1+a13(2)*xin2;
   o13=1/(1+exp(-net13));
   %o13=net13;
   
   net21=a21(1)*o11+a21(2)*o12+a21(3)*o13;	% output layer
   o21=1/(1+exp(-net21));
   %o21=net21;
   net22=a22(1)*o11+a22(2)*o12+a22(3)*o13;
   o22=1/(1+exp(-net22));
   %o22=net22;
   
   e1=1-o21;	% training set 1
   e2=0-o22;
   
   fp=o21*(1-o21);
   a21(1)=a21(1)+mu*e1*fp*o11;	% update output layer
   a21(2)=a21(2)+mu*e1*fp*o12;   
   a21(3)=a21(3)+mu*e1*fp*o13;
   
   fp=o22*(1-o22);
   a22(1)=a22(1)+mu*e2*fp*o11;
   a22(2)=a22(2)+mu*e2*fp*o12;
   a22(3)=a22(3)+mu*e2*fp*o13;
   
   fp=o11*(1-o11);
   d1=fp*(a21(1)*e1 + a22(1)*e2);		% back propagate error
   fp=o12*(1-o12);
   d2=fp*(a21(2)*e1 + a22(2)*e2);
   fp=o13*(1-o13);
   d3=fp*(a21(3)*e1 + a22(3)*e2);
   
   fp=o11*(1-o11);
   a11(1)=a11(1)+mu*d1*fp*xin1;		% hidden layer update
   a11(2)=a11(2)+mu*d1*fp*xin2;
   fp=o12*(1-o12);
   a12(1)=a12(1)+mu*d2*fp*xin1;   
   a12(2)=a12(2)+mu*d2*fp*xin2;   
   fp=o13*(1-o13);
   a13(1)=a13(1)+mu*d3*fp*xin1;   
   a13(2)=a13(2)+mu*d3*fp*xin2;   
   
   % repeat for other class
   
   %xin1=datap2(1,n)/15;
   %xin2=datap2(2,n)/15;
   xin1=m12;
   xin2=m22;
   net11=a11(1)*xin1+a11(2)*xin2;	% hidden layer
   o11=1/(1+exp(-net11));
   %o11=net11;
   net12=a12(1)*xin1+a12(2)*xin2;
   o12=1/(1+exp(-net12));
   %o12=net12;
   net13=a13(1)*xin1+a13(2)*xin2;
   o13=1/(1+exp(-net13));
   %o13=net13;
   
   net21=a21(1)*o11+a21(2)*o12+a21(3)*o13;	% output layer
   o21=1/(1+exp(-net21));
   %o21=net21;
   net22=a22(1)*o11+a22(2)*o12+a22(3)*o13;
   o22=1/(1+exp(-net22));
   %o22=net22;
   
   e1=0-o21;	% training set 2
   e2=1-o22;
   
   fp=o11*(1-o11);
   d1=fp*(a21(1)*e1 + a22(1)*e2);		% back propagate error
   fp=o12*(1-o12);
   d2=fp*(a21(2)*e1 + a22(2)*e2);
   fp=o13*(1-o13);
   d3=fp*(a21(3)*e1 + a22(3)*e2);
   
   fp=o21*(1-o21);
   a21(1)=a21(1)+mu*e1*fp*o11;	% update output layer
   a21(2)=a21(2)+mu*e1*fp*o12;   
   a21(3)=a21(3)+mu*e1*fp*o13;
   
   fp=o22*(1-o22);
   a22(1)=a22(1)+mu*e2*fp*o11;
   a22(2)=a22(2)+mu*e2*fp*o12;
   a22(3)=a22(3)+mu*e2*fp*o13;
   
   fp=o11*(1-o11);
   a11(1)=a11(1)+mu*d1*fp*xin1;		% hidden layer update
   a11(2)=a11(2)+mu*d1*fp*xin2;
   fp=o12*(1-o12);
   a12(1)=a12(1)+mu*d2*fp*xin1;   
   a12(2)=a12(2)+mu*d2*fp*xin2;   
   fp=o13*(1-o13);
   a13(1)=a13(1)+mu*d3*fp*xin1;   
   a13(2)=a13(2)+mu*d3*fp*xin2;   
   
  end; 
end;
   
   
% nets been trained - lets look at output

[X Y] = meshgrid(0:.1:15,0:.1:15);
[Ny Nx]=size(X);
nncls1=zeros(size(X));
nncls2=zeros(size(X));

for nx=1:Nx,
   for ny=1:Ny,
      
      xin1 = X(ny,nx);
      xin2 = Y(ny,nx);
      
      net11=a11(1)*xin1+a11(2)*xin2;	% hidden layer
	   o11=1/(1+exp(-net11));
   	net12=a12(1)*xin1+a12(2)*xin2;
	   o12=1/(1+exp(-net12));
   	net13=a13(1)*xin1+a13(2)*xin2;
	   o13=1/(1+exp(-net13));
   
   	net21=a21(1)*o11+a21(2)*o12+a21(3)*o13;	% output layer
	   o21=1/(1+exp(-net21));
   	net22=a22(1)*o11+a22(2)*o12+a22(3)*o13;
	   o22=1/(1+exp(-net22));
      
      nncls1(ny,nx)=o21;
      nncls2(ny,nx)=o22;
   end;
 end;
 
 figure(2);
 pcolor(X,Y,nncls1);
 shading interp;
 colorbar('vert');
 
 figure(3);
 pcolor(X,Y,nncls2);
 shading interp;
 colorbar('vert');
 