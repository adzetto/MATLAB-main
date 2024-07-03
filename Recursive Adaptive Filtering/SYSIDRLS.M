% sysid program for LMS

fs=1024;
npts=2000;
nstep=8;

np=2;
nz=2;

j=sqrt(-1);

% generate ARMA polynomials

mp1=0.98;
fp1=200;
p1=mp1*cos(2*pi*fp1/fs)+j*mp1*sin(2*pi*fp1/fs);
p1c=mp1*cos(2*pi*fp1/fs)-j*mp1*sin(2*pi*fp1/fs);

A=poly([p1 p1c]);

gain=1;
mz1=0.95;
fz1=400;
z1=mz1*cos(2*pi*fz1/fs)+j*mz1*sin(2*pi*fz1/fs);
z1c=mz1*cos(2*pi*fz1/fs)-j*mz1*sin(2*pi*fz1/fs);

B=gain.*poly([z1 z1c]);

An=zeros(npts/nstep,np+1);	% time arrays
Bn=zeros(npts/nstep,nz+1);
for k=1:npts/nstep,
	An(k,:)=A;
	Bn(k,:)=B;
end;

% generate output data from ZMG input

x=randn([1 npts]);
y=zeros([1 npts]);

for n=1:npts,
	sumB=0;
	for nb=1:nz+1,
		if n-nb > 0 
			sumB=sumB+x(n-nb+1)*B(nb);
		end;
	end;

	sumA=0;
	for na=2:np+1,	% no 0th term...
		if n-na > 0
			sumA=sumA+y(n-na+1)*A(na);
		end;
	end;
    y(n) = sumB - sumA;
end;

% check responses
fftsize=1024;
if fftsize > npts
	fftsize = npts;
end;
Txy=tfestimate(x,y,rectwin(fftsize/2),fftsize/2-1,fftsize);
freq=fs.*(1:fftsize/2)./fftsize;

% compare to polys

Bf=fft(B,fftsize);
Af=fft(A,fftsize);
Hf=Bf./Af;

% rls sysid

nwts=nz+1+np;

% recursive least squares

Id = zeros([nwts nwts]);
P = zeros([nwts nwts]);
H = zeros([nwts 1]);
K = zeros([nwts 1]);
for n=1:nwts,
	Id(n,n)=1;
end;
%Nbar = 1/murel;
Nbar=50;
alpha = (Nbar-1)/Nbar;
beta=1/Nbar;
alphai = 1/ alpha;
pow=std(y)^2;

Ahn=zeros(npts/nstep,np+1);	% time arrays
Bhn=zeros(npts/nstep,nz+1);

% initialization
P = (2.*pow.*Id)^(-1);

% nwts better be 5 for the basis fcn to be right

kstep=1;
for n=nwts+1:npts,
    phi = [x(n) x(n-1) x(n-2) -y(n-1) -y(n-2)];
    K = P*phi'./(alpha + phi*P*phi');
    P = alphai.*(Id - K*phi)*P;
    yp = phi*H;
    ep = y(n) - yp;
    H = H + K*ep;
	kstep=kstep+1;
    if kstep == nstep
		Ahn(round(n/nstep),1)=1;
		Ahn(round(n/nstep),2:3)=H(4:5).';
		Bhn(round(n/nstep),:)=H(1:3).';
		kstep=1;
	end;

end;
Ah=zeros([1 np+1]);
Bh=zeros([1 nz+1]);
Ah(1)=1;
Ah(2:3)=H(4:5).';
Bh(1:3)=H(1:3).';

% compare to RLS

Bfh=fft(Bh,fftsize);
Afh=fft(Ah,fftsize);
Hfh=Bfh./Afh;

% bootstrapped RLS ARMA ID

Nbar=200;
alpha = (Nbar-1)/Nbar;
beta=1/Nbar;
alphai = 1/ alpha;

Abn=zeros(npts/nstep,np+1);	% time arrays
Bbn=zeros(npts/nstep,nz+1);

% initialization
P = (2.*pow.*Id)^(-1);

% nwts better be 5 for the basis fcn to be right
eb=zeros([1 npts]);
ep=zeros([1 npts]);

pow_eb = 0;
pow_ep = 0;
kstep=1;
for n=nwts+1:npts,

    phib = [0 eb(n-1) eb(n-2) -y(n-1) -y(n-2)];
    ypb = phib*H;
    eb(n) = y(n) - ypb;
	pow_eb=alpha*pow_eb + beta*eb(n)^2;

%  if n <= Nbar,
%    alpha = (n-1)/n;
%    alphai = 1/alpha;             % this gives fast startup convergence
%  end;


    phi = [eb(n) eb(n-1) eb(n-2) -y(n-1) -y(n-2)];
    K = P*phi'./(alpha + phi*P*phi');
    P = alphai.*(Id - K*phi)*P;
    yp = phi*H;
    ep(n) = y(n) - yp;
	pow_ep=alpha*pow_ep + beta*ep(n)^2;
    H = H + K*ep(n);

	kstep=kstep+1;
    if kstep == nstep
		Abn(round(n/nstep),1)=1;
		Abn(round(n/nstep),2:3)=H(4:5).';
		Bbn(round(n/nstep),:)=sqrt(pow_eb).*H(1:3).';
		kstep=1;
	end;

end;
Ab=zeros([1 np+1]);
Bb=zeros([1 nz+1]);
Ab(1)=1;
Ab(2:3)=H(4:5).';
Bb(1:3)=sqrt(pow_eb).*H(1:3);

% compare to RLS

Bfb=fft(Bb,fftsize);
Afb=fft(Ab,fftsize);
Hfb=Bfb./Afb;


figure(1);

subplot(2,1,1), plot(freq,(180/pi).*unwrap(angle(Txy(1:fftsize/2))),'-.k');
hold on;
subplot(2,1,1), plot(freq,(180/pi).*angle(Hf(1:fftsize/2)),'--k');
subplot(2,1,1), plot(freq,(180/pi).*angle(Hfh(1:fftsize/2)),'k');
subplot(2,1,1), plot(freq,(180/pi).*angle(Hfb(1:fftsize/2)),':k');
hold off;
ylabel('deg');
title('Frequency Responses');
axis([0 500 -200 200]);

subplot(2,1,2), plot(freq,20.*log10(abs(Txy(1:fftsize/2))),'-.k');
hold on;
subplot(2,1,2), plot(freq,20.*log10(abs(Hf(1:fftsize/2))),'--k');
subplot(2,1,2), plot(freq,20.*log10(abs(Hfh(1:fftsize/2))),'k');
subplot(2,1,2), plot(freq,20.*log10(abs(Hfb(1:fftsize/2))),':k');
hold off;
ylabel('dB');
xlabel('Hz');
legend('FFT Response','Polynomial','RLS System ID','ARMA Bootstrap');
axis([0 500 -30 70]);

itr=nstep.*(1:npts/nstep);
figure(2);
subplot(2,1,1), plot(itr,Bn,'--k',itr,Bhn,'k',itr,Bbn,':k');
title('ARMA Coefficients');
ylabel('B[z]');
subplot(2,1,2), plot(itr,An,'--k',itr,Ahn,'k',itr,Abn,':k');
ylabel('A[z]');
xlabel('Iteration');

