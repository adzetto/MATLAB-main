% Jim Stover's Blend functions


b=.9;
d=.1;

a=1.5;
c=4.5;

dx=.1;
npts=60;
Blnd=zeros([1 npts]);
x=zeros([1 npts]);
for n=1:npts,
   x(n)=n*dx;
   if x(n) <= a
      Blnd(n)=b;
   end;
   if x(n) >= c
      Blnd(n)=d;
   end;
   if x(n) > a
      if x(n) < c
         Blnd(n)=.5*(d+b+(d-b)*sin(pi*((x(n)-a)/(c-a) - .5)));
      end;
   end;
end;

figure(1);
plot(x,Blnd,'k');
title('Information X Blend Function');
ylabel('Confidence in X');
xlabel('Property Value');


% confusion example

b=2;
d=0;

a=-1;
c=1;

dx=.01;
npts=100;
Cnfu=zeros([1 npts]);
x=zeros([1 npts]);
for n=1:npts,
   x(n)=n*dx;
   if x(n) <= a
      Cnfu(n)=b;
   end;
   if x(n) >= c
      Cnfu(n)=d;
   end;
   if x(n) > a
      if x(n) < c
         Cnfu(n)=.5*(d+b+(d-b)*sin(pi*((x(n)-a)/(c-a) - .5)));
      end;
   end;
end;
Cnfu=1 - Cnfu;

figure(2);
plot(x,Cnfu,'k');
title('Confusion Blend Function');
ylabel('Confusion');
xlabel('Property Confidence Difference');

% band-pass blend fcn - body fat mass example

b1=.05;
d1=.95;

b2=d1;
d2=.10

a1=10;
c1=30;

a2=40;
c2=65;

dx=1;
npts=80;
Blndf=zeros([1 npts]);
x=zeros([1 npts]);
for n=1:npts,
   x(n)=n*dx;
   if x(n) <= a1
      Blndf(n)=b1;
   end;
   if x(n) > a1
      if x(n) < c1
         Blndf(n)=.5*(d1+b1+(d1-b1)*sin(pi*((x(n)-a1)/(c1-a1) - .5)));
      end;
   end;
   if x(n) >= c1
      if x(n) < a2
         Blndf(n)=d1;
      end;
   end;
   if x(n) >= a2
      if x(n) < c2
         Blndf(n)=.5*(d2+b2+(d2-b2)*sin(pi*((x(n)-a2)/(c2-a2) - .5)));
      end;
   end;
   if x(n) >= c2
      Blndf(n)=d2;
   end;
end;

figure(3);
plot(x,Blndf,'k');
title('Body Fat Class Blend Example');
ylabel('Confidence');
xlabel('Percent Fat by Weight');
