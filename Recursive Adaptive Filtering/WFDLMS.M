% whitening performance model

fs=1024;  % sample rate

npts=512; % num of samples

% make input-output data

%nzamp = .001;   % use this one to examine convergence of the error spectra
nzamp=10;       % use this one to show model accuracy
nn=1:npts;

num_sigs=3;
f=[80 250 390];
amp=[10 30 20];
x = nzamp.*randn([1 npts]);
y = zeros([1 npts]);

for k=1:num_sigs,
	x = x + amp(k).*(sin(2.*pi.*f(k).*nn./fs));
end;

% okay, now for the filter and output signal

nzeros=4;
zmag=[.99 1.4 .9 .95];
zang=pi.*[.2 .35 .6 .75];
zer=zeros([1 2*nzeros]);
i=sqrt(-1);
for k=1:nzeros,
	zer(2*k-1)=zmag(k)*(cos(zang(k))+i*sin(zang(k)));
	zer(2*k)=zmag(k)*(cos(zang(k))-i*sin(zang(k)));
end; 
Asys=real(poly(zer));
nwts_act = 1 + 2*nzeros;

for n=1:npts,
	y(n)=0;
	for k=1:nwts_act,
		if n-k+1 > 1
			y(n)=y(n)+x(n-k+1)*Asys(k);
		end;
	end;
end;


el = zeros([1 npts]);
er = zeros([1 npts]);
ef = zeros([1 npts]);
eb = zeros([1 npts]);

nwts=2*nwts_act;

% simple LMS Wiener filter

a=zeros([1 nwts]);
pow=std(x)^2;
mumax = 1/(nwts*pow);
murel = 1/(2*nwts);                  % safty margin 

% memory window is approx 1/murel must be at least nwts long!
mu = mumax*murel;

for n=1:npts,
	yp=0;
	for k=1:nwts,
		if n-k+1 > 0
			yp = yp + a(k)*x(n-k+1);
		end;
	end;
	el(n)=y(n)-yp;	

	for k=1:nwts,
		if n-k+1 > 0
			a(k) = a(k) + 2*mu*x(n-k+1)*el(n);    % adaptation
		end;
	end;
end;

% recursive least squares

Id = eye(nwts);
P = zeros([nwts nwts]);
H = zeros([nwts 1]);
H(1)=1;
K = zeros([nwts 1]);
Nbar = 1/murel;
alpha = (Nbar-1)/Nbar;
alphai = 1/ alpha;

% initialization
P = (2.*pow.*Id)^(-1);
phi=zeros([1 nwts]);

for n=1:npts,
  for k=1:nwts,
	if n-k+1 > 0
      phi(k) = x(n-k+1);
    end;
  end;
  K = P*phi'./(alpha + phi*P*phi');
  P = alphai.*(Id - K*phi)*P;
  yp = phi*H;
  er(n)=y(n)-yp;
  H = H + K*er(n);
end;

% freq domain LMS

nfwts=2*nwts;
Aw=zeros([1 nfwts]);

nffts=fix(npts/nfwts);
Px=zeros([1 nfwts]);
for k=1:nffts;
	Px=Px + abs(fft(x((k-1)*nfwts+1:k*nfwts),nfwts)).^2;
end;
Px = Px./nffts;

binave=1;
Pxave=zeros(size(Px));
for k=1:nfwts,
	for j=k-binave:k+binave,
		if j > 0 && j <= nfwts
			Pxave(k) = Pxave(k) + Px(j);
		end;
	end;
end;
mumaxbins = (2*binave+1)./Pxave;

murel = 1;
mubins = mumaxbins.*murel;
af = zeros([1 nfwts]);


overlap=1;
for n=1:npts,
	yp=0;
	for k=1:nfwts/2,
		if n-k+1 > 0
		  yp = yp + af(k)*x(n-k+1);    % prediction
		end;
	end;
	ef(n) = y(n) - yp;

	if overlap > nfwts && n > nfwts
	    ebuf=[zeros([1 nfwts/2]) ef(n-nfwts/2+1:n)];
    	xbuf=x(n-nfwts+1:n);
		Ew=fft(ebuf);
		Xw=fft(xbuf);
		xSpec=conj(Xw).*Ew;
		Aw=Aw + mubins.*xSpec;
		af=real(ifft(Aw));	
		overlap = 0;

        aSpec=conj(Xw).*Xw;
		EoX=xSpec./aSpec;
    
    end;
	overlap = overlap + 1;
end;


figure(1);
subplot(3,1,1), plot(nn,er,'k');
ylabel('RLS Error');
subplot(3,1,2), plot(nn,el,'k');
ylabel('LMS Error');
subplot(3,1,3), plot(nn,ef,'k');
ylabel('FDLMS Error');
xlabel('Iteration');

figure(2);
subplot(3,1,1), plot(20.*log10(Px(1:nwts)),'k');
ylabel('Ref Power');
subplot(3,1,2), plot(20.*log10(Pxave(1:nwts)),'k');
ylabel('Ref Pwr Ave');
subplot(3,1,3), plot(mubins(1:nwts),'k');
ylabel('mu[f]');
xlabel('FFT Bin');



nbins=32;
nstep=16; 
df=.5*fs/nbins;
nlines=fix( (npts-2*nbins)/nstep );
Z=zeros([round(nlines/nstep) nbins]);

for k=1:nlines,
	buf=er( (k-1)*nstep+1:(k-1)*nstep+2*nbins );
	Fbuf=20.*log10(abs(fft(buf))./nbins);
	Z(k,:)=Fbuf(1:nbins);
end;

[X Y] = meshgrid(df:df:.5*fs, 2*nbins+1:nstep:npts);
figure(3);
colormap gray;
mesh(X,Y,Z);
caxis([0 100]);
view(-135,30);
xlabel('Hz');
ylabel('Iteration');
zlabel('dB Error');
title('RLS Error Spectrum');
axis([df .5*fs 1 npts-2*nbins -100 25]);

for k=1:nlines,
	buf=el( (k-1)*nstep+1:(k-1)*nstep+2*nbins );
	Fbuf=20.*log10(abs(fft(buf))./nbins);
	Z(k,:)=Fbuf(1:nbins);
end;

figure(4);
colormap gray;
mesh(X,Y,Z);
caxis([0 100]);
view(-135,30);
xlabel('Hz');
ylabel('Iteration');
zlabel('dB Error');
title('LMS Error Spectrum');
axis([df .5*fs 1 npts-2*nbins -100 25]);

for k=1:nlines,
	buf=ef( (k-1)*nstep+1:(k-1)*nstep+2*nbins );
	Fbuf=20.*log10(abs(fft(buf))./nbins);
	Z(k,:)=Fbuf(1:nbins);
end;

figure(5);
colormap gray;
mesh(X,Y,Z);
caxis([0 100]);
view(-135,30);
xlabel('Hz');
ylabel('Iteration');
zlabel('dB Error');
title('FDLMS Error Spectrum');
axis([df .5*fs 1 npts-2*nbins -100 25]);

figure(6);
zpres=512;
Hw=fft(Asys,zpres);
Hh=fft(af(1:nwts),zpres);
EoXt=ifft(EoX);
EX=fft(EoXt,zpres);
Hep=abs(Hh) + abs(EX);
Hem=abs(Hh) - abs(EX);

ff=fs/zpres:fs/zpres:fs/2;
plot(ff,abs(Hw(1:zpres/2)),'k',ff,abs(Hh(1:zpres/2)),'k--',ff,Hep(1:zpres/2),'k:',ff,Hem(1:zpres/2),'k:');
xlabel('Hz');
ylabel('Transfer Function');
legend('Actual System','Model Response','Error Bounds');

