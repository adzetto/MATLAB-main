% recursive least squares model - WHITENING FILTER

fs=1024;  % sample rate
npts=500; % num of samples - change at npts/2
nhalf=npts/2;

% make input-output data
nn=1:npts;
f0=50;   % 50 and 200 excellent examples

sigamp = 1;
nzamp = .01;
y = nzamp.*randn([1 npts]) + sigamp.*sin(2.*pi.*f0.*nn./fs); 


% simple LMS whitening filter

a=zeros([1 npts]);  % plot array of middle coeff
elms=zeros([1 npts]);
a1p = 0;
a2p = 0;

% model order 3 times x variance -> set mumax
mumax = 1/(2*std(y)^2);
murel = .2;                   % safty margin
mu = mumax*murel;

for n=3:npts,
     yp = -a1p*y(n-1) - a2p*y(n-2);                    % prediction
     ep = y(n) - yp;                                        % pred error
     elms(n)=ep;
     a1p = a1p - 2*mu*y(n-1)*ep;
     a2p = a2p - 2*mu*y(n-2)*ep;
     a(n) = a1p;                                            % save coeff
end;


% recursive least squares whitening filter

h = zeros([1 npts]);
erls = zeros([1 npts]);
Id = zeros([2 2]);
P = zeros([2 2]);
H = zeros([2 1]);
K = zeros([2 1]);
Id(1,1)=1;
Id(2,2)=1;
Nbar = 1/murel;
alpha = (Nbar-1)/Nbar;
alphai = 1/ alpha;

% initialization
phi = [y(2) y(1)];
ypow=std(y)^2;
P = (Nbar.*phi'*phi + .00001.*Id)^(-1);    % slow init
%P=(ypow.*Id)^(-1);                          % fast init

for n=4:npts,
  phi = [y(n-1) y(n-2)];
  K = P*phi'./(alpha + phi*P*phi');
  P = alphai.*(Id - K*phi)*P;
  yp = -phi*H;
  erls(n)=y(n)-yp;
  H = H - K*erls(n);
  h(n)=H(1);
end;

% least squares lattice whitening filter

nwts=2;
Mlat=nwts+2;

d=zeros([Mlat 1]);			% cross cor and variance
rr=zeros([Mlat 1]);

r=zeros([Mlat 1]);			% error sigs
rl=zeros([Mlat 1]);
e=zeros([Mlat 1]);

gmc=zeros([Mlat 1]);		%likelihood and PARCORS
ke=zeros([Mlat 1]);
kr=zeros([Mlat 1]);

A=zeros([Mlat 1]);			% linear pred coeffs
Al=zeros([Mlat 1]);
Bl=zeros([Mlat 1]);

hlat = zeros([1 npts]);
elat = zeros([1 npts]);

gmc(1)=1;					% initialization
ypow=std(y)^2;
%ypow=.00001;
rrl=ypow.*ones([Mlat 1]);
re=ypow.*ones([Mlat 1]);

for n=1:npts,
  e(2)=y(n);
  r(2)=y(n);
  re(2)=alpha*re(2)+y(n)^2;
  rr(2)=re(2);
  for kk=3:Mlat,
    d(kk)=alpha*d(kk)+e(kk-1)*rl(kk-1)/gmc(kk-2);
    ke(kk)=d(kk)/re(kk-1);
	kr(kk)=d(kk)/rrl(kk-1);
	e(kk)=e(kk-1)-kr(kk)*rl(kk-1);
	r(kk)=rl(kk-1)-ke(kk)*e(kk-1);
	re(kk)=re(kk-1)-d(kk)*kr(kk);
	rr(kk)=rrl(kk-1)-d(kk)*ke(kk);
	gmc(kk-1)=gmc(kk-2)-(rl(kk-1)^2)/rrl(kk-1);
  end;

% Levinson recursion

  Al(3)=-kr(3);
  Bl(3)=-ke(3);
  if nwts > 1
	for ii=4:Mlat,
	  Al(ii)=-kr(ii);
	  Bl(ii)=-ke(ii);
	  pp=ii-1;
	  for kk=3:pp,
		A(kk)=Al(kk)-kr(ii)*Bl(ii-kk+2);
		Bl(ii-kk+2)=Bl(ii-kk+2)-ke(ii)*Al(kk);
		Al(kk)=A(kk);
	  end
	end
	for ii=1:nwts,
	  A(ii)=Al(ii+2);	% shift em back down by to avoid confusion
	end;
  end

% update time-lagged variables

  for kk=2:Mlat,
    rrl(kk)=rr(kk);
	rl(kk)=r(kk);
  end;
  hlat(n)=A(1);
  elat(n)=e(Mlat);
end;


figure(1);
plot(nn,a,'k-.',nn,h,'k',nn,hlat,'k:');
xlabel('Iteration');
ylabel('a1 Coefficient');
legend('LMS','RLS','Lattice');
axis([0 50 -2 2]);

figure(2);
plot(nn,elms,'k-.',nn,erls,'k',nn,elat,'k:');
xlabel('Iteration');
ylabel('error signal');
legend('LMS','RLS','Lattice');
axis([0 50 -1 1]);

X = [ (nzamp^2+.5*sigamp^2) (sigamp^2*cos(2*pi*f0/fs)); (sigamp^2*cos(2*pi*f0/fs)) (nzamp^2+.5*sigamp^2)];
[V,D]=eig(X);
