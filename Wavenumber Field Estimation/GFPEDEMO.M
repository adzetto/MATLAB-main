% GFPE by Dave Swanson, May 29, 1997
% benchmark agreement confirmed...finally!
% k-space integration for surface wave
% ref: Gilbert & Di, JASA 94(4), Oct 1993

j=sqrt(-1);
t_start=cputime;

f = 100;	% frequency in Hz
dr=20;		% range step
rmax=2000;	% max down range
zs=1;		% source height  
zr=1;		% receiver height
c0=344.;		% base sound speed
lambda=c0/f;

cturblev=0.01;
%cturblev=0.0;

npts=1024;	% set npts and zmax
zmax=1024;	% height in m

zattn=300;
nattn=fix(npts*zattn/zmax);

%Zg = 38.79 + j*38.41; 	% @ 10 Hz
Zg = 12.81 + j*11.62;   % 4-parm @ 100Hz, Attenborough, JSV 99
%Zg = 5.96+j*2.46;		% @ 1kHz

% dr/dz > about 80 for good results

dz=zmax/npts;

drdz=dr/dz
lambdz=lambda/dz

dk=2*pi/npts;

Nrech=fix(zr/dz)+1;	% n offset of reciever height

Nr=round(rmax/dr);	%number of range steps
r=dr:dr:Nr*dr;

z=0:dz:zmax-dz;
kp=0:dk:2*pi-dk;
kp=kp./dz;		% needed for Nyquist scaled physical wavenumber

kb=zeros(size(kp));
oz=ones(size(z));	% z is npts long - same as kp

cz=zeros(size(z));
%cz=c0.*oz + z./10;		% ken's downward profile
cz=c0.*oz - z./10;		% ken's upward profile
%cz=c0.*oz;				% constant profile
%for n=1:npts,			% down-up-constant profile
%	if z(n) <= 100
%		cz(n)=c0+z(n)/10;
%	elseif z(n) <= 300
%		cz(n)=c0+10-(z(n)-100)/10;
%	elseif z(n) > 300
%		cz(n)=cz(n-1);
%	end;
%end;

kr = 2*pi*f/c0;			% kr is constant with range
beta = kr/Zg;			% surface wavenumber
beta_np2=beta*ones([1 npts/2]);

kb(1:npts/2)=kp(1:npts/2) - beta_np2;
for i=npts:-1:npts/2+1,
	kb(i)=(i-npts-1)*dk/dz - beta; % -k space
end;

kz=(2*pi*f)./cz;

% Gaussian starting field at r = 0

psi0=sqrt(kr/2).*(exp(-((.5*kr).*(z-zs.*oz)).^2) + ...
	((Zg-1)/(Zg+1)).*exp(-((.5*kr).*(z+zs.*oz)).^2));

% Since Zg, c0, and kr are constant with range we can compute
% a few vectors only once

Rkp=zeros([1 npts]);
Rkpm=zeros([1 npts]);
dk2=zeros([1 npts]);

Rkp=(kp-beta.*oz)./(kp+beta.*oz);	% reflection factors
Rkpm=(-kp-beta.*oz)./(-kp+beta.*oz); % -k ref factor

%dk2=kz.^2-(kr^2).*oz;			% wave number variation
%e_vari=exp(j.*dr.*dk2./(2*kr));

% k-delta sq now includes attn layer
e_vari=zeros([1 npts]);
alpha=zeros([1 npts]);
top=dz*npts/2;
for n=1:npts/2,
	if n > nattn
		alpha(n)=(z(n)-zattn)/(top-zattn);
		alpha(n)=alpha(n)^2;
	end;
	k1=kz(n)+j*alpha(n);
	Carg=j*dr*(k1*k1-kr*kr)/(2*kr);	% dk2...w/alpha
	e_vari(n)=exp(Carg);
end;
e_vari(npts/2+1:npts)=zeros([1 npts/2]);  % zero pad

%
% Here are the wave number filters
%
e_step=exp( (j*dr).*(sqrt((kr^2).*oz - kp.^2)-kr.*oz) );
e_step_b=exp(j*dr*(sqrt(kr^2 - beta^2)-kr));
e_jbz=exp(-j.*beta.*z);

% now fill a npts high by Nr wide grid

psi=zeros([1 npts]);
psif=zeros([1 npts]);
u=zeros([1 npts]);
psi_temp=zeros([1 npts]);
P_grid=zeros([nattn Nr]);

t_init=cputime;

psi=psi0;		% starting field at r=0

cturb=zeros(size(e_vari));

for nr=1:Nr,

	psi(1)=.5*psi(1);		% trap integration
	psif=fft(psi);			% FFT

	T=0;
	for i=1:npts/4,			% k-space integration
		T = psif(i)/kb(i) + psif(npts-i+1)/kb(npts-i+1) +T;
	end;
	B=(-dk*2*beta*T*e_step_b/(2*pi)).*e_jbz;  % ifft 

	u=zeros(size(psif));
	for n=2:npts/2,			% Xio Shuffle
		t1=psif(n);
		t2=psif(npts+2-n);
		u(n)=(t1 + t2*Rkp(n)).*e_step(n);
		u(npts+2-n)=(t2 + t1*Rkpm(n)).*e_step(n);
	end;
	u(1)=psif(1)*(1 + Rkp(1)).*e_step(1);
	u(npts/2+1)=psif(npts/2+1)*(1 + Rkp(npts/2+1)).*e_step(npts/2+1);

	A=ifft(u);				% i FFT	

	cturb=exp(j.*pi.*cturblev.*randn(size(e_vari)));								
	psi_temp=cturb.*e_vari.*(A + B);	% direct & reflected summed

	P_grid(1:nattn,nr)=psi_temp(1:nattn).'./sqrt(r(nr));

	psi=psi_temp;					% get fft's ready for next loop

end;

[X Y] = meshgrid(r,z(1:nattn));

figure(1);

pcolor(X,Y, 20.*log10(abs(P_grid(1:nattn,:))));
colormap gray;
shading interp;
ylabel('Height m','color','k');
xlabel('Range m','color','k');
title_bufdB=sprintf('dB f=%.f Zs=%.1f Zr=%.1f dz=%.4f dr=%.2f N=%d c0=%.f',f,zs,zr,dz,dr,npts,c0);
title(title_bufdB,'color','k');

