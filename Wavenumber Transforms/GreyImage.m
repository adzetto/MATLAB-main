% open grayscale jpeg image
% this is just a plain old jpeg image on disk - any jpeg will work
IM = imread('House Dec 2004 gs.jpg');

info = imfinfo('House Dec 2004 gs.jpg');  % more info on pic if you want it

% get the number of rows and columns, find center, and capture a 256 x 256
% image to work with (makes FFT processing easy)
[nr, nc]=size(IM);
cnr=round(nr/2);
cnc=round(nc/2);
gi=zeros([256,256]);
gi=IM(cnr-128:cnr+127,cnc-128:cnc+127);

% check on subimage
figure(1);
imagesc(gi);
colormap(gray);
imagesc(gi);        % I use imagesc here to automatically scale brightness
title('Original Image');

% this does a 2DFFT and shifts the <0,0> wavenumber to the center
gifft=fftshift( fft2(gi) );

% display shifted 2DFFT
figure(2);
colormap(gray);
imagesc(20.*log10(abs(gifft)));
title('FFT of Original Image');

% this is a crude low-pass filter (easy to calculate)
han=hanning(256);
% note that raising to the 8th power increases the emphasis
lpwfilt=(han.^8)*(han.^8)';

hpwfilt=ones(size(lpwfilt))-lpwfilt;
hpwfilt=hpwfilt.^8; % this just enhances the hp effect even more

% simple point by point low pass filter in wavenumber domain!
lpgifft=gifft.*lpwfilt;

% this is the resulting image - un-shift first, then ifft, the real part
% the imaginary part can have residual junk in it and we don't need it
lpgi=real(ifft2(fftshift(lpgifft)));

% show the low-pass filtered image
figure(3);
colormap(gray);
imagesc(lpgi);
title('Low Pass Filtered Image');

% repeat for the high pass filtered image
hpgifft=gifft.*hpwfilt;
hpgi=real(ifft2(fftshift(hpgifft)));
figure(4);
colormap(gray);
imagesc(hpgi);
title('High Pass Filtered Image');

% show the low pass FFT - log scale to show dynamics
figure(5);
colormap(gray);
imagesc(20.*log10(abs(lpgifft)+1));
title('Low Pass Image FFT');

% here's an unsharp operation just for fun - easy in frequency domain
sharpgi=hpgi+.3.*lpgi;
figure(6);
colormap(gray);
imagesc(sharpgi);
title('HP + .3 LP Image');

% show high pass FFT
figure(7);
colormap(gray);
imagesc(20.*log10(abs(hpgifft)+1));
title('High Pass Image FFT');

% OK here's the "just for fun" stuff
sl=8;
funfft=gifft;
% zero out the bands
funfft(128-sl:128+sl,:)=zeros(2*sl+1,256);
funfft(:,128-sl:128+sl)=zeros(256,2*sl+1);
figure(8);
colormap(gray);
imagesc(20.*log10(abs(funfft)+1));
title('Just for Fun FFT');

% show the underlaying image
funim=real(ifft2(fftshift(funfft)));
figure(9);
colormap(gray);
imagesc(real(funim));
title('Just fo Fun Image');

% now for some more fun
sl=8;
fun2fft=zeros(256,256);
fun2fft(128-sl:128+sl,1:128-sl-1)=gifft(128-sl:128+sl,1:128-sl-1);
fun2fft(128-sl:128+sl,128+sl+1:256)=gifft(128-sl:128+sl,128+sl+1:256);
fun2fft(:,128-sl:128+sl)=gifft(:,128-sl:128+sl);

figure(10);
colormap(gray);
imagesc(20.*log10(abs(fun2fft)+1));
title('More Fun FFT');

fun2im=real(ifft2(fftshift(fun2fft)));
figure(11);
colormap(gray);
imagesc(real(fun2im));
title('More Fun Image');

