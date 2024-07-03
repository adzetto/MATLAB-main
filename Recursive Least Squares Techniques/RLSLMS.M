% recursive least squares model

fs=1024;  % sample rate
f1=100;        % starting zero
f2=256;        
f3=400;        % ending zero

npts=500; % num of samples - change at npts/2

% make input-output data
nn=1:npts;
f=25;          % 25 amp 5 works great
               % 256 shows no effect
%f=256;

sigamp = 5;
nzamp = 1;

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
     yp = x(n) + a1p*x(n-1) + a2p*x(n-2);   % prediction
     ep = y(n) - yp;                        % pred error
     a0p = a0p + 2*mu*x(n)*ep;              % update a0   
     a1p = a1p + 2*mu*x(n-1)*ep;            % update a1
     a2p = a2p + 2*mu*x(n-2)*ep;            % update a2
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
%  if n <= Nbar,
%    alpha = (n-1)/n;
%    alphai = 1/alpha;             % this gives fast startup convergence
%  end;
  phi = [x(n) x(n-1) x(n-2)];
  K = P*phi'./(alpha + phi*P*phi');
  P = alphai.*(Id - K*phi)*P;
  yp = phi*H;
  ep = y(n) - yp;
  H = H + K*ep;
  hrls(n) = H(2);
  erls(n) = ep;
end;

figure(1);
plot(nn,hlms,'k-.',nn,hrls,'k',nn,htrue,'k:');
legend(' LMS',' RLS', 'True');
xlabel('Iteration');
ylabel('h1 Coefficient');
axis([0 500 -2.5 +2.5]);

figure(2);
plot(nn,elms,'k-.',nn,erls,'k');
legend(' LMS',' RLS');
xlabel('Iteration');
ylabel('error signal');
axis([0 500 -10 10]);

X = [ (nzamp^2+.5*sigamp^2) (sigamp^2*cos(2*pi*f/fs)); (sigamp^2*cos(2*pi*f/fs)) (nzamp^2+.5*sigamp^2)];
[V,D]=eig(X);
