% adaptive null former for linear array
% Eigenvalue solution

nels=8;
nsnaps=500;

T=.1;	%time overlap for FFT's
f=100;
c=343;
lambda=c/f;
d=1.25;
w=2*pi*f;
kd=w*d/c;

nsrc=3;
Ampsrc=[2 3 1];
Angsrc=[45 90 135];	%degrees

nsrc=1;        % comment for 3 sources otherwise 1 src
Ampsrc=3;
Angsrc=60;	%degrees

nzl=1e-4;	% e-7 smallest, e+1 largest
nzsc=1;

srcang=Angsrc*pi/180.;

Y=zeros([nsnaps nels]);
A=zeros([nels 1]);

j=sqrt(-1);

for n=1:nsnaps,
	phases=nzsc.*(2*pi).*([rand rand rand]);  % scr's random at each snapshot
	for m=1:nels,
		Y(n,m)=0;
		for i=1:nsrc,
			Y(n,m)=Y(n,m)+Ampsrc(i)*exp(j*(w*T*n+phases(i)))*exp(-j.*kd.*(m-1).*cos(srcang(i)));
		end;
		Y(n,m)=Y(n,m) + nzl*(randn + j*randn);
	end;
end;

% auto correlation matrix - gotta love Matlab!
R=Y'*Y;

[V, D] = eig(R);    % eigenvectors and eigenvalues

% check order of eigenvalues - signals last in new versions of Matlab
[Dmax Imax]=max(D);

music=zeros(size(V(:,1)));
for m=1:nels-nsrc,
	music=music + V(:,m);   % this should be the noise only eigenvectors
end;
music=music./(nels-nsrc);

npts=1024;
musicF=fftshift(fft(music,npts));

% all eigenvectors in freq domain
EF=zeros(npts,nels);
for m=1:nels,
	EF(:,m)=fftshift(fft(V(:,m),npts));
end;

piscale=-pi:2*pi/npts:pi*(npts-1)/npts;
thetascale=(180/pi).*real(acos(piscale./kd));

% plot the loudest source (last eigenvector)
figure(1);
plot(thetascale,20.*log10(abs(EF(:,nels))),'k');
if nsrc == 3
    title('Signal Eigenvector Responses');
    hold on;
    plot(thetascale,20.*log10(abs(EF(:,nels-1))),'k--');
    plot(thetascale,20.*log10(abs(EF(:,nels-2))),'k:');
    hold off;
    l1=sprintf('\\lambda_1=%.3f',D(nels,nels)/(nsnaps*nels));
    l2=sprintf('\\lambda_2=%.3f',D(nels-1,nels-1)/(nsnaps*nels));
    l3=sprintf('\\lambda_3=%.3f',D(nels-2,nels-2)/(nsnaps*nels));
    legend(l1,l2,l3);
else
    buff=sprintf('Eigenvalue is %.3f',D(nels,nels)./(nsnaps*nels));
    title(buff);
end;
ylabel('Beam Response dB');
xlabel('degrees');
h=gca;
set(h,'XTick',[0 45 90 135 180]);
axis([0 180 -120 20]);

% plot all the other eigenvectors
figure(2);
plot(thetascale,20.*log10(abs(EF(:,1))),'k:');
hold on;
for m=2:nels-nsrc,
	plot(thetascale,20.*log10(abs(EF(:,m))),'k:');
end;
hold off;
ylabel('Beam Response dB');
xlabel('degrees');
title('Noise Eigenvectors');
h=gca;
set(h,'XTick',[0 45 90 135 180]);
axis([0 180 -120 20]);

% plot the sum of the noise eigenvectors (MUSIC)
figure(3);
plot(thetascale,20.*log10(abs(musicF)),'k');
ylabel('Beam Response dB');
xlabel('degrees');
title('MUSIC Beam Response');
h=gca;
set(h,'XTick',[0 45 90 135 180]);
axis([0 180 -120 20]);

% steer the beam and hold the nulls
steerangle=input('Steering angle?');
steertheta=pi*(steerangle)/180;
Steer=zeros([nels 1]);
for m=1:nels,
	Steer(m)=exp(j*kd*(m-1)*cos(steertheta));
end;

% this is what all the fuss is about!
Rinv=inv(R);
denom=Steer'*Rinv*Steer;
gmv=(Rinv*Steer)./denom;
gmvs=Steer;

AFp=fftshift(fft(gmv,npts));
SFp=fftshift(fft(gmvs,npts));

figure(4);
plot(thetascale,20.*log10(abs(AFp)),'k');
hold on;
plot(thetascale,20.*log10(abs(SFp)),'k:');
hold off;
legend('MVDR','Bartlett');
ylabel('Beam Response');
xlabel('degrees');
buff=sprintf('MVDR and Bartlet Beams Steered to %f. Degrees',steerangle);
title(buff);
h=gca;
set(h,'XTick',[0 45 90 135 180]);
axis([0 180 -120 20]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lets optimize our beams to the sources
% psuedo reduced-order method

Zmusic=roots(music);
ZV=Zmusic;
ZVm=abs(ZV);
ZVa=angle(ZV);
epsilon=.005;	% the detection threshold
i=sqrt(-1);
srcdet=ones(size(ZV));
for m=1:nels-1,				% check zeros
	if ZVm(m) > 1+epsilon
		ZVm(m) = 10.;		% push out if > 1
		srcdet(m)=0;
	end;
	if ZVm(m) < 1-epsilon
		ZVm(m) = .1;		% push in if < 1
		srcdet(m)=0;
	end;
	ZV(m)=ZVm(m)*(cos(ZVa(m))+i*sin(ZVa(m)));
end;

mymusic=poly(ZV);
scale=0;
for m=1:nels,
	scale=scale+abs(mymusic(m))^2;
end;
scale=sqrt(scale)*2;
mymusic=mymusic./scale;
mymusicF=fftshift(fft(mymusic,npts));

figure(5);
plot(thetascale,20.*log10(abs(mymusicF)),'k');
ylabel('Beam Response');
xlabel('degrees');
title('Pseudo Reduced-Order MUSIC Beam Response');
h=gca;
set(h,'XTick',[0 45 90 135 180]);
axis([0 180 -120 20]);

% how many sources?
srcZ=zeros([nels-1,nels-1]);
srcbeams=zeros([nels,nels-1]);
estsrcang=zeros([1 nels-1]);
srcoffset=zeros([1 nels-1]);
numdetsrc=0;
for m=1:nels-1,
	if srcdet(m)==1			% source here
		numdetsrc=numdetsrc+1;
		srcZ(:,numdetsrc)=ZV;
		srcZ(m,numdetsrc)=0.;	% push source zero to 0
		estsrcang(numdetsrc)=(180/pi)*acos(angle(ZV(m))/kd);
		srcoffset(numdetsrc)=m;
	end;
end;
numdetsrc	% show the user
estsrcang   % show user src angs

lastfig=5;

if numdetsrc > 0
	figure(lastfig);
	for m = 1:numdetsrc,
		ZSV=srcZ(:,m);
		srcbeams(:,m)=poly(ZSV).';  % make a beam
		scale=0;
		for mm=1:nels,
			scale=scale+abs(srcbeams(mm,m))^2;
		end;
		scale=sqrt(scale)*2;
		srcbeams(:,m)=srcbeams(:,m)./scale;
		srcbeamF=fftshift(fft(srcbeams(:,m).',npts));
		if m == 1
			plot(thetascale,20.*log10(abs(srcbeamF)),'k');
			ylabel('Beam Response');
			xlabel('degrees');
			title('Pseudo Reduced-Order Signal Beam Response');
			h=gca;
			set(h,'XTick',[0 45 90 135 180]);
			axis([0 180 -120 20]);
			hold on;
		end;
		if m == 2
			plot(thetascale,20.*log10(abs(srcbeamF)),'k:');
		end;
		if m == 3
			plot(thetascale,20.*log10(abs(srcbeamF)),'k--');
		end;
	end;
	hold off;
end;

% ok get serious! when looking at src x
% we want zeros at the other sources
% and everywhere else!

nz=nels-numdetsrc+1;	% # of zeros to work with + 1
myang=2*pi/nz;		% spreading angle
myZVang=zeros([1 nz]);
for m=1:nz,
	myZVang(m)=m*myang-myang*nz/2;
end;

myZVang
myZVang=acos(myZVang/kd)	% pre shift them

myZV=zeros([nels-1 nels-1]);
i=sqrt(-1);
if numdetsrc > 0
	for m=1:numdetsrc,
	
		% remove the one closest to our source
		% shuffle the rest
	    closest=-1;
		maxdiffz=999;
		for mm=1:nz-1,
			diffz=abs(myZVang(mm)-estsrcang(m)*(pi/180));
			if diffz < maxdiffz
				closest=mm;
				maxdiffz=diffz;
			end;
		end;
		for mm=closest:nz-1,	% shuffle
			myZVang(mm)=myZVang(mm+1);
		end;			
		% assemble zeros
		for mm=1:nz-2,
			ang=kd*cos(myZVang(mm));
			myZV(mm,m)=cos(ang)+i*sin(ang);
		end;
		% add other src zeros
		offset=nz-1;
		for mm=1:numdetsrc,
			if mm ~= m
				ang=kd*cos(estsrcang(mm)*(pi/180));
				myZV(offset,m)=cos(ang)+i*sin(ang);
				offset=offset+1;
			end;
		end;
	end;

	% okay now we have all the zeros!

	lastfig=lastfig+1;
	figure(lastfig);

	for my = 1:numdetsrc,
		ZSV=myZV(:,my);
		srcbeams(:,my)=poly(ZSV).';  % make a beam
		scale=0;
		for mm=1:nels,
			scale=scale+abs(srcbeams(mm,my))^2;
		end;
		scale=sqrt(scale)*2;
		srcbeams(:,my)=srcbeams(:,my)./scale;
		srcbeamF=fftshift(fft(srcbeams(:,my).',npts));
		if my == 1
			plot(thetascale,20.*log10(abs(srcbeamF)),'k');
			ylabel('Beam Response');
			xlabel('degrees');
			title('Null Synthesis Signal Beam Response');
			h=gca;
			set(h,'XTick',[0 45 90 135 180]);
			axis([0 180 -120 20]);
			hold on;
		end;
		if my == 2
			plot(thetascale,20.*log10(abs(srcbeamF)),'k:');
		end;
		if my == 3
			plot(thetascale,20.*log10(abs(srcbeamF)),'k--');
		end;
	end;
	hold off;
end;
