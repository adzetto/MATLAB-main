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

gain=5;
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


% LMS sys id

e=zeros([1 npts]);

nah=2;
nbh=2;
Ah=0.01.*randn([1 nah+1]);
Bh=0.01.*randn([1 nbh+1]);
Ah(1)=1;
Ahn=zeros(npts/nstep,np+1);	% time arrays
Bhn=zeros(npts/nstep,nz+1);

pow_x=std(x)^2;
pow_y=std(y)^2;
mu_B=.1/(nbh*pow_x);
mu_A=.01/(nah*pow_y);

kstep=1;
for n=1:npts,
	sumBh=0;
	for k=1:nbh+1,
		if n-k >= 0 
			sumBh=sumBh+x(n-k+1)*Bh(k);
		end;
	end;

	sumAh=0;
	for i=2:nah+1,	% no 0th term...
		if n-i >= 0
			sumAh=sumAh+y(n-i+1)*Ah(i);
		end;
	end;
    yh = sumBh - sumAh;
    e(n) = y(n) - yh;

% LMS updates

	for k=1:nbh+1,
		if n-k >= 0 
			Bh(k)=Bh(k)+2*mu_B.*x(n-k+1)*e(n);
		end;
	end;

	for i=2:nah+1,	% no 0th term...
		if n-i >= 0
			Ah(i)=Ah(i)-2*mu_A.*y(n-i+1)*e(n);
		end;
	end;
	kstep=kstep+1;
    if kstep == nstep
		Ahn(round(n/nstep),:)=Ah;
		Bhn(round(n/nstep),:)=Bh;
		kstep=1;
	end;
end;

% compare to LMS

Bfh=fft(Bh,fftsize);
Afh=fft(Ah,fftsize);
Hfh=Bfh./Afh;

% Let's see if we can bootstrap is (ARMA signal model)

eb=zeros([1 npts]);
ebe=zeros([1 npts]);

nab=2;
nbb=2;
Ab=0.01.*randn([1 nab+1]);
Bb=0.01.*randn([1 nbb+1]);
Ab(1)=1;
Abn=zeros(npts/nstep,np+1);	% time arrays
Bbn=zeros(npts/nstep,nz+1);

mu_Ab=.0125/(nah*pow_y);
mu_Bb=.01/(nah*pow_y);;

pow_eb=0;
Nbar=500;
alpha=(Nbar-1)/Nbar;
beta=1/Nbar;

kstep=1;
for n=1:npts,

% bootstrap
	sumBh=0;
	for k=2:nbb+1,
		if n-k >= 0 
			sumBh=sumBh+eb(n-k+1)*Bb(k);
		end;
	end;

	sumAh=0;
	for i=2:nab+1,	% no 0th term...
		if n-i >= 0
			sumAh=sumAh+y(n-i+1)*Ab(i);
		end;
	end;
    yh = sumBh - sumAh;
    eb(n) = y(n) - yh;
	sumBh=sumBh+eb(n)*Bb(1);	% updated after bootstrap
    yh = sumBh - sumAh;
    ebe(n) = y(n) - yh;			% bootstrap error

	pow_eb=alpha*pow_eb + beta*eb(n)^2;

% LMS updates

	for k=1:nbh+1,
		if n-k > 0 
			Bb(k)=Bb(k)+2*mu_Bb*eb(n-k+1)*ebe(n);
		end;
	end;

	for i=2:nah+1,	% no 0th term...
		if n-i > 0
			Ab(i)=Ab(i)-2*mu_Ab*y(n-i+1)*ebe(n);
		end;
	end;
	kstep=kstep+1;
    if kstep == nstep
		Abn(round(n/nstep),:)=Ab;
		Bbn(round(n/nstep),:)=sqrt(pow_eb).*Bb;
		kstep=1;
	end;
end;
Bb=sqrt(pow_eb).*Bb;


% compare to LMS

Bfb=fft(Bb,fftsize);
Afb=fft(Ab,fftsize);
Hfb=Bfb./Afb;

itr=nstep.*(1:npts/nstep);
figure(2);
subplot(2,1,1), plot(itr,Bn,'--k',itr,Bhn,'k',itr,Bbn,':k');
title('ARMA Coefficients');
ylabel('B[z]');
h=subplot(2,1,2), plot(itr,An,'--k',itr,Ahn,'k',itr,Abn,':k');
ylabel('A[z]');
xlabel('Iteration');
legend(h,'Polynomial','LMS System ID','ARMA Bootstrap');


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
legend('FFT Response','Polynomial','LMS System ID','ARMA Bootstrap');
axis([0 500 -30 70]);

% this is a tricky plot!
% An&Bn are the true coeffs, Ahn&Bhn are the LMS est, Abn&Bbn are the
% LMS bootstrapped  all 3 coeffs are plotted in each set
itr=nstep.*(1:npts/nstep);
figure(2);
subplot(2,1,1), plot(itr,Bn,'--k',itr,Bhn,'k',itr,Bbn,':k');
title('ARMA Coefficients');
ylabel('B[z]');
subplot(2,1,2), plot(itr,An,'--k',itr,Ahn,'k',itr,Abn,':k');
ylabel('A[z]');
xlabel('Iteration');
axis([0 npts -1.5 1.5]);

