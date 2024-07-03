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


Nbar=200;
alpha = (Nbar-1)/Nbar;

% lattice sysid only for order 2 systems

% initialization - double direct lattice

u = [y; x];
Rrue0 = (1/npts).*u*u';
Rrur0 = Rrue0;

Kur1=zeros([2 2]);
Kue1=zeros([2 2]);
rup1l=zeros([2 1]);
ru1l=zeros([2 1]);
Rrue1 = Rrue0;
Rrur1 = Rrue0;

Kur2=zeros([2 2]);
Kue2=zeros([2 2]);
rup2l=zeros([2 1]);
ru2l=zeros([2 1]);
Rrue2 = Rrue0;
Rrur2 = Rrue0;

z21 = [0; 0];
rup0 = z21;
rup1 = z21;
rup2 = z21;
ru0 = z21;
ru1 = z21;
ru2 = z21;
I2 = [1 0; 0 1];

Ahn=zeros(npts/nstep,np+1);	% time arrays
Bhn=zeros(npts/nstep,nz+1);

kstep=1;
for n=1:npts,
	rupl0 = rup0;						% time shift
	rupl1 = rup1;
	rupl2 = rup2;
	rul0 = ru0;
	rul1 = ru1;
	rul2 = ru2;

	eu0 = u(:,n);							% stage 0
    ru0 = eu0;
	eup0 = eu0;
	rup0 = ru0;
    Rrue0 = alpha.*Rrue0 + eu0*eu0';
	Rrur0 = alpha.*Rrur0 + ru0*ru0';

	eup1 = eup0 - Kur1*rupl0;			% stage 1
	rup1 = rupl0 - Kue1*eup0;
	Kue1 = Kue1 + rup1*eu0'/Rrue0;
	Kur1 = Kur1 + eup1*rul0'/Rrur0;

	eu1 = eu0 - Kur1*rul0;
	ru1 = rul0 - Kue1*eu0;
	Rrue1 = alpha.*Rrue1 + eu1*eup1';
	Rrur1 = alpha.*Rrur1 + ru1*rup1';

	eup2 = eup1 - Kur2*rupl1;			% stage 2
	rup2 = rupl1 - Kue2*eup1;
	Kue2 = Kue2 + rup2*eu1'/Rrue1;
	Kur2 = Kur2 + eup2*rul1'/Rrur1;

	eu2 = eu1 - Kur2*rul1;
	ru2 = rul1 - Kue2*eu1;
	Rrue2 = alpha.*Rrue2 + eu2*eup2';
	Rrur2 = alpha.*Rrur2 + ru2*rup2';

% Levinson recursion

	A2 = -Kur2;
	A1 = -Kur1 + Kur2*Kue1;
	A0 = [1 -sqrt((Rrue2(1,1)/Rrue2(2,2))); ...
  		  -sqrt((Rrue2(2,2)/Rrue2(1,1))) 1];

	kstep = kstep+1;
	if kstep == nstep
		Ahn(round(n/nstep),:) = [A0(1,1) A1(1,1) A2(1,1)];
		Bhn(round(n/nstep),:) = [-A0(1,2) -A1(1,2) -A2(1,2)];	
		kstep = 1;
	end;
end;

% compare to polys

Bfh=fft([-A0(1,2) -A1(1,2) -A2(1,2)],fftsize);
Afh=fft([A0(1,1) A1(1,1) A2(1,1)],fftsize);
Hfh=Bfh./Afh;

% bootstrapped lattice


Nbar=200;				% 200 work great!
alpha = (Nbar-1)/Nbar;

% lattice sysid only for order 2 systems

% initialization - double direct lattice

pow = std(y)^2;;
Rrue0 = [pow 0; 0 pow];
Rrur0 = Rrue0;

Kur1=zeros([2 2]);
Kue1=zeros([2 2]);
rup1l=zeros([2 1]);
ru1l=zeros([2 1]);
Rrue1 = Rrue0;
Rrur1 = Rrue0;

Kur2=zeros([2 2]);
Kue2=zeros([2 2]);
rup2l=zeros([2 1]);
ru2l=zeros([2 1]);
Rrue2 = Rrue0;
Rrur2 = Rrue0;

z21 = [0; 0];
rup0 = z21;
rup1 = z21;
rup2 = z21;
ru0 = z21;
ru1 = z21;
ru2 = z21;
I2 = [1 0; 0 1];

Abn=zeros(npts/nstep,np+1);	% time arrays
Bbn=zeros(npts/nstep,nz+1);

kstep=1;
for n=1:npts,
	rupl0 = rup0;						% time shift
	rupl1 = rup1;
	rupl2 = rup2;
	rul0 = ru0;
	rul1 = ru1;
	rul2 = ru2;

	eub0 = [y(n); 0];
	eub1 = eub0 - Kur1*rul0;
	eub2 = eub1 - Kur2*rul1;

	eu0 = [y(n); eub2(1)];					% stage 0 bootstrapped
    ru0 = eu0;
	eup0 = eu0;
	rup0 = ru0;
    Rrue0 = alpha.*Rrue0 + eu0*eu0';
	Rrur0 = alpha.*Rrur0 + ru0*ru0';

	eup1 = eup0 - Kur1*rupl0;			% stage 1
	rup1 = rupl0 - Kue1*eup0;
	Kue1 = Kue1 + rup1*eu0'/Rrue0;
	Kur1 = Kur1 + eup1*rul0'/Rrur0;

	eu1 = eu0 - Kur1*rul0;
	ru1 = rul0 - Kue1*eu0;
	Rrue1 = alpha.*Rrue1 + eu1*eup1';
	Rrur1 = alpha.*Rrur1 + ru1*rup1';

	eup2 = eup1 - Kur2*rupl1;			% stage 2
	rup2 = rupl1 - Kue2*eup1;
	Kue2 = Kue2 + rup2*eu1'/Rrue1;
	Kur2 = Kur2 + eup2*rul1'/Rrur1;

	eu2 = eu1 - Kur2*rul1;
	ru2 = rul1 - Kue2*eu1;
	Rrue2 = alpha.*Rrue2 + eu2*eup2';
	Rrur2 = alpha.*Rrur2 + ru2*rup2';

% Levinson recursion

	A2 = -Kur2;
	A1 = -Kur1 + Kur2*Kue1;
	A0 = [1 -1; -1 1];

	kstep = kstep+1;
	if kstep == nstep
		Abn(round(n/nstep),:) = [A0(1,1) A1(1,1) A2(1,1)];
		Bbn(round(n/nstep),:) = (sqrt(Rrue2(1,1)/Nbar)).*[-A0(1,2) -A1(1,2) -A2(1,2)];	
		kstep = 1;
	end;
end;

% compare to polys

Bfb=fft((sqrt(Rrue2(1,1)/Nbar)).*[-A0(1,2) -A1(1,2) -A2(1,2)],fftsize);
Afb=fft([A0(1,1) A1(1,1) A2(1,1)],fftsize);
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
legend('FFT Response','Polynomial','RLS Lattice System ID','ARMA Lattice Bootstrap');
axis([0 500 -30 70]);

itr=nstep.*(1:npts/nstep);
figure(2);
subplot(2,1,1), plot(itr,Bn,'--k',itr,Bhn,'k',itr,Bbn,':k');
title('ARMA Coefficients');
ylabel('B[z]');
subplot(2,1,2), plot(itr,An,'--k',itr,Ahn,'k',itr,Abn,':k');
ylabel('A[z]');
xlabel('Iteration');

