% Helmholz-Huygens geometry

xs=-50;
ys=-50;
xr=50;
yr=-50;

x=[-300:1:300];
[nr nc]=size(x);

nz=0;
skip=nc/20;
nskips=0;
yval=nz*randn;
for n=1:nc,
  y(n)=yval;
  nskips=nskips+1;
  if nskips >= skip
	yval = nz*randn;
    nskips=0;
  end;
end;

xi=zeros(size(x));
yi=zeros(size(x));


for n=1:nc,
  theta=atan2(xr-x(n),y(n)-yr);
  rs=sqrt((x(n)-xs)^2 + (y(n)-ys)^2);
  xi(n)=x(n)-rs*sin(theta);
  yi(n)=y(n)+rs*cos(theta);
end;



figure(1);
plot(x,y,'y',xs,ys,'wo',xr,yr,'wx',xi,yi,'g');
axis equal;
