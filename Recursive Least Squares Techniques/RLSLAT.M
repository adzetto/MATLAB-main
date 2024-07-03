% recursive least squares model

% Dave Swanson, Dec 8, 1996
% runs LMS, RLS and Lattice Wiener filters
% filter zero moves in steps for tracking analysis

fs=1024;  % sample rate
f1=100;        % starting zero
f2=256;        
f3=400;        % ending zero

npts=500; % num of samples - change at npts/2

% make input-output data
nn=1:npts;
f=25;          % f = 25 amp 5 works great
               % f = 256 shows no effect

sigamp = 5;
nzamp = 1;

% input is white noise plus a sinusoid (to widen eigenvalues)

x = nzamp.*randn([1 npts]) + sigamp.*(sin(2.*pi.*f.*nn./fs)); 
y = zeros([npts 1]);
a1 = -2*cos(2*pi*f1/fs);
a2 = -2*cos(2*pi*f2/fs);
a3 = -2*cos(2*pi*f3/fs);

htrue=zeros([1 npts]);
hlms=zeros([1 npts]);    % plot array of middle coeff
hrls = zeros([1 npts]);
elms=zeros([1 npts]);
erls=zeros([1 npts]);

y(1)=x(1);
y(2)=x(2)+a1*x(1);
htrue(1)=a1;
htrue(2)=a1;
n3=fix(npts/3);
for n=3:n3,
 y(n) = x(n) + a1*x(n-1) + x(n-2);
 htrue(n)=a1;
end;
for n=n3+1:npts-n3,
 y(n) = x(n) + a2*x(n-1) + x(n-2);
 htrue(n)=a2;
end;
for n=npts-n3+1:npts,
 y(n) = x(n) + a3*x(n-1) + x(n-2);
 htrue(n)=a3;
end;

% simple LMS sys id

a0p = 0;
a1p = 0;
a2p = 0;

% model order 3 times x variance -> set mumax
pow=std(x)^2;
mumax = 1/(2*pow);
murel = .05;                  % safty margin
mu = mumax*murel;

for n=3:npts,
     yp = x(n) + a1p*x(n-1) + a2p*x(n-2);    % prediction
     ep = y(n) - yp;                         % pred error
     a0p = a0p + 2*mu*x(n)*ep;
     a1p = a1p + 2*mu*x(n-1)*ep;
     a2p = a2p + 2*mu*x(n-2)*ep;
     hlms(n) = a1p;                                    % save coeff
    elms(n) = ep;                                      % save error
end;

% recursive least squares

Id = zeros([3 3]);
P = zeros([3 3]);
H = zeros([3 1]);
K = zeros([3 1]);
Id(1,1)=1;
Id(2,2)=1;
Id(3,3)=1;
Nbar = 1/murel;
alpha = (Nbar-1)/Nbar;
alphai = 1/ alpha;

% initialization
P = (2.*pow.*Id)^(-1);

for n=3:npts,
  if n <= Nbar,
    alpha = (n-1)/n;
    alphai = 1/alpha;             % this gives fast startup convergence
  end;
  phi = [x(n) x(n-1) x(n-2)];
  K = P*phi'./(alpha + phi*P*phi');
  P = alphai.*(Id - K*phi)*P;
  yp = phi*H;
  ep = y(n) - yp;
  H = H + K*ep;
  hrls(n) = H(2);
  erls(n) = ep;
end;


% least squares lattice whitening filter

nwts=2;
Mlat=nwts+2;

d=zeros([Mlat 1]);			% cross cor and variance
dg=zeros([Mlat 1]);			% cross cor and variance
rr=zeros([Mlat 1]);

r=zeros([Mlat 1]);			% error sigs
rp=zeros([Mlat 1]);			% a priori
rl=zeros([Mlat 1]);
rlp=zeros([Mlat 1]);		% a priori
e=zeros([Mlat 1]);
ep=zeros([Mlat 1]);			% a priori
eg=zeros([Mlat 1]);
egp=zeros([Mlat 1]);		% a priori

gmc=zeros([Mlat 1]);		%likelihood and PARCORS
ke=zeros([Mlat 1]);
kr=zeros([Mlat 1]);
kg=zeros([Mlat 1]);

A=zeros([Mlat 1]);			% linear pred coeffs
Al=zeros([Mlat 1]);
Bl=zeros([Mlat 1]);

hlat = zeros([1 npts]);
elat = zeros([1 npts]);

rr=pow.*ones([Mlat 1]);
rrl=pow.*ones([Mlat 1]);
re=pow.*ones([Mlat 1]);

gmc(1)=1;

Lb=zeros([(nwts+1) (nwts+1)]);
for ii=1:nwts+1,
	Lb(ii,ii)=1;
end;
hvect=zeros([(nwts+1) 1]);
gvect=zeros([(nwts+1) 1]);

for n=1:npts,
  e(2)=x(n);
  r(2)=x(n);
  ep(2)=e(2);
  rp(2)=r(2);

  re(2)=alpha*re(2)+x(n)^2;
  rr(2)=re(2);

  egp(2)=y(n)-kg(2)*rp(2);

  dg(2)=alpha*dg(2)+y(n)*r(2);
  kg(2)=dg(2)/rr(2);

  eg(2)=y(n)-kg(2)*r(2);

  for kk=3:Mlat,

% the a priori error updates needed to demo numerical error (next 2)
	 ep(kk)=ep(kk-1)-kr(kk)*rlp(kk-1);
	 rp(kk)=rlp(kk-1)-ke(kk)*ep(kk-1);

    d(kk)=alpha*d(kk) + rl(kk-1)*e(kk-1)/gmc(kk-2);
    ke(kk)=d(kk)/re(kk-1);
	kr(kk)=d(kk)/rrl(kk-1);

	e(kk)=e(kk-1)-kr(kk)*rl(kk-1);
	r(kk)=rl(kk-1)-ke(kk)*e(kk-1);

	re(kk)=re(kk-1)-d(kk)*kr(kk);
	rr(kk)=rrl(kk-1)-d(kk)*ke(kk);
	gmc(kk-1)=gmc(kk-2)-(rl(kk-1)^2)/rrl(kk-1);


% switch kg and you'll have numerical error!
	 dg(kk)=alpha*dg(kk)+eg(kk-1)*r(kk)/gmc(kk-1); % a bit noisy!
	 kg(kk)=dg(kk)/rr(kk);

%	 egp(kk)=egp(kk-1)-kg(kk)*rp(kk);				% real clean!
%     kg(kk)=kg(kk)+egp(kk)*r(kk)/rr(kk);			

	eg(kk)=eg(kk-1)-kg(kk)*r(kk);

  end;

% Levinson recursion

  Al(3)=-kr(3);
  Bl(3)=-ke(3);
  Lb(1,2)=Bl(3);
  if nwts > 1
	for ii=4:Mlat,
	  Al(ii)=-kr(ii);
	  Bl(ii)=-ke(ii);
	  Lb(ii-3,ii-1)=Bl(ii);
	  pp=ii-1;
	  for kk=3:pp,
		A(kk)=Al(kk)-kr(ii)*Bl(ii-kk+2);
		Bl(ii-kk+2)=Bl(ii-kk+2)-ke(ii)*Al(kk);
		Al(kk)=A(kk);
  	    Lb(ii-kk+1,ii-1)=Bl(ii-kk+2);
	  end
	end
    gvect(1)=kg(2);
	for ii=1:nwts,
	  A(ii)=Al(ii+2);	% shift em back down by 2 to avoid confusion
      gvect(ii+1)=kg(ii+2);
	end;
  end

% update time-lagged variables

  for kk=2:Mlat,
    rrl(kk)=rr(kk);
	rl(kk)=r(kk);
    rlp(kk)=rp(kk);
  end;
  hvect=Lb*gvect;		% Lb is really L^T
  hlat(n)=hvect(2);
  elat(n)=eg(Mlat);
end;




figure(1);
plot(nn,hlms,'k:',nn,hrls,'k',nn,hlat,'k-.',nn,htrue,'k--');
legend(' LMS',' RLS','Lattice', 'True');
xlabel('Iteration');
ylabel('h1 Coefficient');
axis([0 500 -2.5 +2.5]);

figure(2);
plot(nn,elms,'k:',nn,erls,'k',nn,elat,'k-.');
legend(' LMS',' RLS','Lattice');
xlabel('Iteration');
ylabel('error signal');
axis([0 500 -10 10]);

X = [ (nzamp^2+.5*sigamp^2) (sigamp^2*cos(2*pi*f/fs)); (sigamp^2*cos(2*pi*f/fs)) (nzamp^2+.5*sigamp^2)];
[V,D]=eig(X);
