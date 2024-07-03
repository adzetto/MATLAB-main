% DAC examples with oversampling reconstruction

fs=400;
f=70;
f8=8*fs;        % this fs is for up to 8x oversampling
T=1/f8;
npts=200;

y=zeros(size([1,npts]));
yd=zeros(size([1,npts]));
ys=zeros(size([1,npts]));
n=[1:npts];

for i=1:npts/8,
	k=8*(i-1)+1;
	y(k)=cos(2*pi*f*T*(k-1));   % analog waveform
	yd(k)=y(k);                 % sample and hold
	ys(k)=y(k);                 % zero padded signal
	for l=1:7,
		k=k+1;
		y(k)=cos(2*pi*f*T*(k-1));
		yd(k)=yd(k-1);
		ys(k)=0;
	end;
end;

figure(1);
plot(n,y,'k');
hold on
plot(n(1:8:npts),yd(1:8:npts),'k.');    %every 8th point
plot(n,yd,'k:');
hold off;
xlabel('Sample Number');

% make a 0th, 1st(2x)3rd(4x) and 7th(8x order oversampling interpolator
h0=zeros(size([1,64]));
h1=zeros(size([1,64]));
h2=zeros(size([1,64]));
h3=zeros(size([1,64]));

for i=1:8,
	h0(i)=1;            % 1x
end;
h1=.125.*conv(h0,h0);   % 2x integral of h0 is 8
h2=.125.*conv(h1,h1);   % 4x
h3=.125.*conv(h2,h2);   % 8x

% plot the interpolation filters
figure(2);
plot(h0,'ko-');
hold on;
plot(h1,'k*-');
plot(h2,'ks-');
plot(h3,'k^-');
axis([0 64 -.25 1.25]);
hold off;
legend('1x','2x','4x','8x');
xlabel('Filter Coefficient');
ylabel('Response');
title('Interpolation filters for 1x, 2x, 4x, and 8x Over Sampling');

% convolve the oversampled signal with each interpolation filter and plot
yr0=conv(ys,h0);
yr1=conv(ys,h1);
yr2=conv(ys,h2)./sqrt(max(h2));
yr3=conv(ys,h3)./sqrt(max(h3));

% shift by conv delay
figure(3);
plot(y(1:150-32),'k');
hold on
plot(yr0(1:150-32),'ko:');
plot(yr1(8:150-32+8),'k*:');    %delay is 8
plot(yr2(15:150-32+15),'ks:');  %delay is 8*2-1
plot(yr3(29:150-32+29),'k^:');  % delay is 15*2-1
axis([0 90 -1 1]);
hold off;
legend('Original signal','1x Recon','2x Recon','4x Recon','8x Recon');
xlabel('Sample Number');
ylabel('Response');
title('Comparison of Over Sampled Signal Reconstruction');

figure(4);
plot(n,ys,'k');
hold on
plot(n(1:8:npts),ys(1:8:npts),'k.');    %every 8th point
hold off;
xlabel('Sample Number');

