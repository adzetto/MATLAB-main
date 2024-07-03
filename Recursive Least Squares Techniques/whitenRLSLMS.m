% recursive least squares model - WHITENING FILTER

fs=1024;  % sample rate
npts=500; % num of samples - change at npts/2
nhalf=npts/2;

% make input-output data
nn=1:npts;
f0=50;   % 50 and 256 excellent examples
%f0=200;

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
P = (Nbar.*phi'*phi + .00001.*Id)^(-1);

for n=4:npts,
  phi = [y(n-1) y(n-2)];
  K = P*phi'./(alpha + phi*P*phi');
  P = alphai.*(Id - K*phi)*P;
  yp = -phi*H;
  erls(n)=y(n)-yp;
  H = H - K*erls(n);
  h(n)=H(1);
end;

figure(1);
plot(nn,a,'k-.',nn,h,'k');
xlabel('Iteration');
ylabel('a1 Coefficient');
legend('LMS','RLS');

figure(2);
plot(nn,elms,'k-.',nn,erls,'k');
xlabel('Iteration');
ylabel('error signal');
legend('LMS','RLS');

X = [ (nzamp^2+.5*sigamp^2) (sigamp^2*cos(2*pi*f0/fs)); (sigamp^2*cos(2*pi*f0/fs)) (nzamp^2+.5*sigamp^2)];
[V,D]=eig(X);
