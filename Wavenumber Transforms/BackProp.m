% back prop simultation
% same as catscan.m except traces summed in catf
load headf;

figure(1);
colormap(gray);
imagesc(20.*log10(abs(headf)+1));

catf=zeros([256 256]);

numscan = 128;
dth=pi/numscan;

% we are summing the 2DFFT traces again for convenience.  In a real
% backprop algorithm, the raw scans are summed, then low pass filtered
for n=1:numscan,
 	th = (n-1)*dth;
	tanth=tan(th);
	tanthc=tan(.5*pi-th);
	if th < pi/4
		for m=1:256,
			r=m;
			c=round(128+(m-128)*tanth);
			catf(r,c)=catf(r,c)+headf(r,c); % note the sum
		end;
	elseif th < (3*pi/4)
		for m=1:256,
			c=m;
			r=round(128+(m-128)*tanthc);
			catf(r,c)=catf(r,c)+headf(r,c);
		end;
	elseif th < pi
		for m=1:256,
			r=m;
			c=round(128-(129-m)*tanth);
			catf(r,c)=catf(r,c)+headf(r,c);
		end;
	end;
end;

% here the low frequency wavenumbers near the center are summed, making
% them brighter than they should be in the image
figure(2);
colormap(gray);
%axis('off');
set(gcf,'color','w');
imagesc(20.*log10(abs(catf)+1));

% the back propagated image is foggy because the low frequencies are added
% more than the high frequencies due to overlap near the origin
headcat=ifft2(fftshift(catf));
figure(3);
colormap(gray);
%axis('off');
set(gcf,'color','w');
imagesc(abs(headcat));

% this is a simple low-pass filter in the frequency domain based on the
% amplitude of the wavenumber, making it zero at the center
[X,Y]=meshgrid(-128:127, -128:127);
hpwfilt=((X.^2 + Y.^2).^.5)./(sqrt(2*128^2));
% high pass filter the 2DFFT back propagation sum
hpcatf = catf.*hpwfilt;
% inverse 2DFFT
hphead=ifft2(fftshift(hpcatf))+ 128.*ones([256 256]);
figure(4);
imagesc(abs(hphead));
colormap(gray);
%axis('off');
set(gcf,'color','w');


