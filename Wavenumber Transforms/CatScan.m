% cat scan simultation

% "headf.mat" is a line drawing of a brain cross section 256 x 256 pixels
% already a 2DFFT and saved as a mat file
load headf;

figure(1);
colormap(gray);
imagesc(20.*log10(abs(headf)+1));
%axis('off');
set(gcf,'color','w');

catf=zeros([256 256]);

numscan = 128;
dth=pi/numscan;
% To "simulate a "scan" we could define each scan line of sight and sum
% across the head image (ifft2 of headf) pixels, but a much easier way is
% to just sample the FFT along the normal to the scan direction.  This is
% significantly faster than rescanning each pixel over and over again for
% each scan direction.  If you don't believe me look a numscan=1


%axis('off');
set(gcf,'color','w');
title(titlebuf);% collect the bins that correspond to the scan angle in catf
titlebuf=sprintf('%d Scans',numscan);
for n=1:numscan,
 	th = (n-1)*dth;% + (40/180)*pi; % the commented term is for a 40 deg single scan
	tanth=tan(th);
	tanthc=tan(.5*pi-th);
	if th < pi/4
		for m=1:256,
			r=m;
			c=round(128+(m-128)*tanth);
			catf(r,c)=headf(r,c);       % here low frequencies can be overwritten
		end;
	elseif th < (3*pi/4)
		for m=1:256,
			c=m;
			r=round(128+(m-128)*tanthc);
			catf(r,c)=headf(r,c);
		end;
	elseif th < pi
		for m=1:256,
			r=m;
			c=round(128-(129-m)*tanth);
			catf(r,c)=headf(r,c);
		end;
	end;
end;

% here we show the 2DFFT of the scan (really just the trace from the 2DFFT
figure(2);
colormap(gray);
imagesc(20.*log10(abs(catf)+1));

% now just do the inverse 2DFFT to see what the scans produced
headcat=ifft2(fftshift(catf));
figure(3);
colormap(gray);
imagesc(abs(headcat));
%axis('off');
set(gcf,'color','w');
title(titlebuf);

