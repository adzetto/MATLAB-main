% maximal length spread spectrum signal generator
% Dave Swanson, 21 August, 1996, Experimental test version for Doug

% n=7, 127 element chip sequence [7,1]s
% just set n as 11,7,6,5,4,3, or 2
n=3;						% n=7 is denoted as [7,1]s
nchips=2.^n - 1;
chips=zeros([1 nchips]);
%rand('seed',0);
%reg=round(rand(1, 11));		% randomize initial registers
reg=[1 0 1 0 1 0 0 0 0 0 0];	% fixed start (same results for comparisons)

for nn=1:nchips,
  chips(nn)=reg(n);			% works for n=2,3,4,5,6,7,11
  reg(11)=reg(10);
  reg(10)=reg(9);
  reg(9)=reg(8);
  reg(8)=reg(7);
  reg(7)=reg(6);
  reg(6)=reg(5);
  reg(5)=reg(4);
  reg(4)=reg(3);
  reg(3)=reg(2);
  reg(2)=reg(1);
  reg(1)=reg(1)+chips(nn);
  if reg(1) > 1,			% modulo 2 sum needed
	 reg(1)=0;
  end;	
end;
% chips is an MLS 0 - 1 sequence
% other convenient sequences - Dixon, Table 3.7
%	n		2^n-1		sequence
%	2		3			[2,1]s
%	3		7			[3,1]s
%	4		15			[4,1]s
%	5		31			[5,1]s
%	6		63			[6,1]s
%	7		127			[7,1]s
%	9		511			[9,4]s
%	10		1023		[10,3]s
%	11		2047		[11,1]s
%	15		32767		[15,1]s


chips2=chips.*2-1;	% make it a +-1 sequence

%Yf=fft(chips2);
%Yfa=Yf.*conj(Yf);
%ya=real(ifft(Yfa));	% check to see if circular autocor is ok
%plot(real(ya(1:20));

% ok lets make a sinewave in a npts vector

npts=16384;
fs=1024;		% sample rate
c=350;			% speed of sound
mlsres=c/fs;	% system resolution in m/bin
%f0=11;
%f0=64;			% carrier freq in Hz
f0=200;
%f0=205;
arg=2*pi*f0/fs;
nt=1:npts;

Amp=30000;
y=Amp.*sin(arg.*nt);

% now set the chip rate and make biphase modulation vector
% the same length as y repeating MLS sequence if necessary

c=350;					% speed of sound
%chiprate=4;			% chips per second
%chiprate=16;			% chips per second
%chiprate=33;			% chips per second
chiprate=100;			% chips per second
chiprateT=1/chiprate;	% chip period in seconds
chiplength=c*chiprateT	% length of 1 chip in m

mlsrate=chiprate/nchips;
mlsrateT=1/mlsrate;
mlslength=c*mlsrateT	% length of mls sequence in m

Ts=1/fs;				% sample period in seconds

% scroll through the chips to make a vector same length as y
time=Ts.*nt;			% time of each sample
biphasemod=zeros(size(y));	% modulation vector
chipoffset=1;			% pointer to current chip in sequence
modoffset=1;			% pointer into modulation vector sample

for k=1:npts,
	if time(k) < modoffset*chiprateT,
		biphasemod(k)=chips2(chipoffset);
	else
		chipoffset=chipoffset+1;
		modoffset=modoffset+1;
		if chipoffset > nchips,
			chipoffset=1;
		end;
		biphasemod(k)=chips2(chipoffset);
	end;
end;

yss=y.*biphasemod;		% spread spectrum signal f0+=chiprate in main lobe
%yss=nz;				% white noise in case you want to compare

nppts=400;
figure(1);
subplot(3,1,1),plot(time(1:nppts),y(1:nppts),'k');
title('Spread Spectrum Time Signal');
buf=sprintf('%.f Hz',f0);
ylabel(buf);

subplot(3,1,2),plot(time(1:nppts),biphasemod(1:nppts),'k');
ylabel('MLS Sequence');
axis([0 time(nppts) -2 2]);

subplot(3,1,3),plot(time(1:nppts),yss(1:nppts),'k');
xlabel('sec');
buf=sprintf('%.f Hz SS',f0);
ylabel(buf);

fftsize=512;			% fft size

fftbuflength=c*Ts*fftsize	% fft buffer length in m
Pyf=spectrum(yss,fftsize);	% Pyf(:,1) is power spec
							% Pyf(:,2) is confidence interval
nave=fix(npts/fftsize);		% number of averages done

% plot the spread spectrum signal out
freqscl=linspace(0,((.5*fftsize-1)/(.5*fftsize))*fs/2,fftsize/2);

figure(2);
plot(freqscl,10.*log10(Pyf(1:fftsize/2,1)),'k');
xlabel('Hz');
ylabel('dB');
title('Transmitted Spreadspectrum');


% create one big FFT to do precise propagation (time delay)

Yf=fft(yss);

% "duct" prop model

Xd1=33.;			% path 1 dist
w=linspace(-pi*fs,+pi*fs*(npts-1)/npts,npts);
jkx1=-j.*w.*Xd1./c;

Xd2=38.;			% path 2 dist 
jkx2=-j.*w.*Xd2./c;

prop=fftshift(exp(jkx1) + exp(jkx2));	% shift to matlab format

% prop spread spectrum data, invert FFT and add noise
Yfprop=Yf.*prop;
Nzlev=0;
yssprop=ifft(Yfprop)+Nzlev.*randn(size(yss));

spec=0;				% set to 1 for Hanning window, 0 rectangular window 
					% 2 for circular correlation corrected
if spec == 1
	P=spectrum(yss, yssprop,fftsize);
	Txy=P(:,4).';	% get transfer fcn
end;

if spec == 0		% I may modify this to correct for circular correlation...
  Txy=zeros([1 fftsize]);
  xspec=zeros([1 fftsize]);
  for nf=1:nave,
	in=yss((1+(nf-1)*fftsize):(nf*fftsize));
	out=yssprop((1+(nf-1)*fftsize):(nf*fftsize));
	inF=fft(in);
	outF=fft(out);
	Txy=Txy+outF./inF;		% straightforward averages transfer fcn
	xspec=xspec+outF.*conj(inF);
  end;
  Txy=Txy./nave;
  xspec=xspec./nave;	
end;

if spec == 2		% to correct for circular correlation...
  Txy=zeros([1 fftsize]);
  for nf=1:nave,
    in=zeros([1 fftsize]);
	in(1+.5*fftsize:fftsize)=yss((1+(nf-.5)*fftsize):((nf)*fftsize));	% half buffer
	out=yssprop((1+(nf-1)*fftsize):(nf*fftsize));
	inF=fft(in);
	outF=fft(out);
	Txy=Txy+outF./inF;		% straightforward averages transfer fcn
  end;
  Txy=Txy./nave;
end;


% ifft to get response vs distance
xdata=abs(ifft(Txy));

% determine distance window to plot response
xscale=linspace(0,fftbuflength*(fftsize-1)/fftsize,fftsize);
xstart=.95*Xd1;
xstop=1.05*Xd2;
nstart=fix(xstart/(c*Ts));
nstop=round(xstop/(c*Ts));
if nstop > .5*fftsize
	nstop=.5*fftsize;
end;

nstart=1;
nstop=fftsize/2;

figure(3);
plot(xscale(nstart:nstop),xdata(nstart:nstop),'k');
xlabel('distance in m');
ylabel('magnitude');
title('Transfer Fcn Distance Measurement');

figure(4);
subplot(2,1,1), plot(freqscl,(180/pi).*angle(Txy(1:fftsize/2)),'k');
axis([0 512 -180 180]);
title('Channel Transfer Fcn');
ylabel('deg');
subplot(2,1,2), plot(freqscl,20.*log10(abs(Txy(1:fftsize/2))),'k');
ylabel('dB');
axis([0 512 -30 30]);

xcorr=abs(ifft(xspec));
figure(5);
plot(xscale(nstart:nstop),xcorr(nstart:nstop),'k');
xlabel('distance in m');
ylabel('magnitude');
title('Correlation Distance Measurement');

