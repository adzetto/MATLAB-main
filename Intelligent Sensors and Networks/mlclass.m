% ML classifier model

npts=500;
m11=7.5;
s11=3;
m21=8;
s21=1;

m12=7;
s12=2;
m22=7;
s22=1;

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

% check results
m1est1=mean(datap1(1,:));
m2est1=mean(datap1(2,:));

%datap1(1,:)=datap1(1,:)-m11;
%datap1(2,:)=datap1(2,:)-m21;
%datap2(1,:)=datap2(1,:)-m12;
%datap2(2,:)=datap2(2,:)-m22;

%covar1=datap1*datap1';	% de-mean em and get the covariance
%covar1=(1/npts).*covar1;
%covar2=datap2*datap2';
%covar2=(1/npts).*covar2;

covar1=R1*[s11^2 0; 0 s21^2]*R1';
covar2=R2*[s12^2 0; 0 s22^2]*R2';

th1=.5*atan2(2*covar1(1,2),(covar1(1,1)-covar1(2,2)));
thdeg1=th1*180/pi

rotback1=[cos(th1) sin(th1); -sin(th1) cos(th1)];
covard1=rotback1*covar1*rotback1';
s1est1=sqrt(covard1(1,1));
s2est1=sqrt(covard1(2,2));

[m11 s11 m21 s21 tilt1*(180/pi)]
[m1est1 s1est1 m2est1 s2est1 thdeg1]


[X Y] = meshgrid(0:.1:15,0:.1:15);

Xp1=X-m11;
Yp1=Y-m21;
Sk1=det(covar1);
P01=1/(2*pi*sqrt(Sk1));
ic1=covar1^(-1);
maha1=0.5.*(Xp1.*ic1(1,1)+Yp1.*ic1(2,1)).*Xp1 + (Xp1.*ic1(1,2)+Yp1.*ic1(2,2)).*Yp1;

Xp2=X-m12;
Yp2=Y-m22;
Sk2=det(covar2);
P02=1/(2*pi*sqrt(Sk2));
ic2=covar2^(-1);
maha2=0.5.*(Xp2.*ic2(1,1)+Yp2.*ic2(2,1)).*Xp2 + (Xp2.*ic2(1,2)+Yp2.*ic2(2,2)).*Yp2;

pdf2d1=P01.*exp(-maha1);
pdf2d2=P02.*exp(-maha2);

figure(2);
pcolor(X,Y,pdf2d1);
shading interp;
xlabel('x_1');
ylabel('x_2');

figure(3);
pcolor(X,Y,pdf2d2);
shading interp;
xlabel('x_1');
ylabel('x_2');

% define sum of training set
ms1=(m11+m12)*.5;
ms2=(m21+m22)*.5;

datap1(1,:)=datap1(1,:)-ms1;
datap1(2,:)=datap1(2,:)-ms2;
datap2(1,:)=datap2(1,:)-ms1;
datap2(2,:)=datap2(2,:)-ms2;
covars=(.5/npts).*(datap1*datap1' + datap2*datap2');
Xps=X-ms1;
Yps=Y-ms2;
Ss=det(covars);
P0s=1/(2*pi*sqrt(Ss));
ics=covars^(-1);
mahas=0.5.*(Xps.*ics(1,1)+Yps.*ics(2,1)).*Xps + (Xps.*ics(1,2)+Yps.*ics(2,2)).*Yps;
pdfsum=P0s.*exp(-mahas);
sumthresh=P0s*.0111;	% 3 sigma bound, 2 sig is .1353

figure(4);
pcolor(X,Y,pdfsum);
shading interp;
xlabel('x_1');
ylabel('x_2');

pdfmax=zeros(size(pdf2d2));
[Ny Nx]=size(X);
pmax=0;
for ny=1:Ny,
   for nx=1:Nx,
      pdfmax(nx,ny)=0;
      if pdf2d1(nx,ny)>=pdf2d2(nx,ny)
         pdfmax(nx,ny)=pdf2d1(nx,ny);
      else
         pdfmax(nx,ny)=pdf2d2(nx,ny);
      end;
      if pdfsum(nx,ny)<=sumthresh
         pdfmax(nx,ny)=0.;
      end;   
      if pdfmax(nx,ny) > pmax
         pmax=pdfmax(nx,ny);
      end;
   end;
end;

figure(5);
surf(X,Y,pdfmax);
shading interp;
xlabel('x_1');
ylabel('x_2');
zlabel('pdf');

% find the discriminant function

edge=abs(del2(pdfmax)).^.125;

figure(6);
pcolor(X,Y,edge);
shading interp;
colormap (bone);
xlabel('x_1');
ylabel('x_2');
title('Edge Detection of Composite Probability');

figure(7);
pstart=log10(0.001*pmax);
pstop=log10(pmax);
levels=logspace(pstart,pstop,10);
contour(X,Y,pdfmax,levels,'k');
xlabel('x_1');
ylabel('x_2');
