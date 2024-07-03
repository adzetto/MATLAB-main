% open grayscale jpeg image

IM = imread('House Dec 2004 gs.jpg');

info = imfinfo('House Dec 2004 gs.jpg');
figure(1);
imshow(IM);
colormap(gray);
[nr, nc]=size(IM);

% boxcar low pass filter
N=64;
Hbc=(1/N^2).*ones([N,N]);
%IMbc=imfilter(IM,Hbc);
IMbc=filter2(Hbc,IM);
figure(2);
%imshow(IMbc)
imagesc(IMbc);
axis off;
colormap(gray);

% x and y gradients
[IMdx, IMdy]=gradient(single(IM));
% add 127 to cover negative numbers
IMdx=IMdx+127.*ones(size(IMdx));  
IMdy=IMdy+127.*ones(size(IMdy));

Hdx=[1 0 -1; 2 0 -2; 1 0 -1];
IMdxf=filter2(Hdx,IM);
Hdy=[-1 -2 -1; 0 0 0; 1 2 1];
IMdyf=filter2(Hdy,IM);

% lets zoom in to see some detail
nx1=round(.2*nr);
nx2=round(.5*nr);
ny1=round(.2*nc);
ny2=round(.5*nc);

figure(3);
imagesc(IMdxf(nx1:nx2,ny1:ny2));
axis off;
colormap(gray);
title('+x-dir derivative');

figure(4);
imagesc(IMdyf(nx1:nx2,ny1:ny2));
axis off;
colormap(gray);
title('+y-dir derivative');

% Sobel edge detection
IMsob=abs(IMdxf)+abs(IMdyf);
maxB=max(max(IMsob));
IMsobinv=maxB.*ones(size(IMsob))-IMsob;    % this just inverts it
figure(5);
imagesc(IMsobinv(nx1:nx2,ny1:ny2));
axis off;
colormap(gray);
title('Inverted Sobel Edges approx by Kirsh Operator');

%Note: the 0.25 scale factor is used to compare the two image without
%autoscaling that is done in imagesc()

% LaPlacian with diagonals added in
Hlap=[-1 -1 -1; -1 8 -1; -1 -1 -1];
IMlap=filter2(Hlap,IM)+127.*ones(size(IM));
figure(6);
image(.25.*IMlap(nx1:nx2,ny1:ny2));
axis off;
colormap(gray);
title('LaPlacian');

% LaPlacian of Gaussian (Matlabs)
HLoG=[-1/6 -2/3 -1/6; -2/3 4+1/3 -2/3; -1/6 -2/3 -1/6];
IMLoG=filter2(HLoG,IM);
figure(7);
image(.25.*IMLoG(nx1:nx2,ny1:ny2));
axis off;
colormap(gray);
title('LaPlacian of Gaussian');

% LaPlacian of Gaussian
sig=.66;
[X,Y]=meshgrid(-1:1,-1:1);
R=X.^2+Y.^2;
tsq=2*sig^2;
pi=acos(-1);
W=(1/(pi*sig^4)).*(ones(size(X))-R./tsq).*exp(-R./tsq);
IMms=filter2(W,IM);
figure(8);
image(.6.*IMms(nx1:nx2,ny1:ny2));
title('LaPlacian of Gaussian \sigma=0.66');
colormap(gray);
axis off;

figure(9);
subplot(1,2,1),image(.25.*IM(nx1:nx2,ny1:ny2));
axis off;
colormap(gray);
title('Original');
subplot(1,2,2),image(.6.*IMms(nx1:nx2,ny1:ny2));
axis off;
colormap(gray);
title('LoG \sigma=0.66');


Hunsharp=fspecial('unsharp');
IMs=imfilter(IM,Hunsharp,'replicate');
figure(10);
imagesc(IMs(nx1:nx2,ny1:ny2));
axis off;
colormap(gray);
title('Matlab Unsharp');

Hmotion=fspecial('motion',25,15);
IMs=imfilter(IM,Hmotion,'replicate');

Hprewitt=fspecial('prewitt');
IMp=imfilter(IM,Hprewitt,'replicate');
IMp2=imfilter(IM,Hprewitt','replicate');

Hsobel=fspecial('sobel');
IMs=imfilter(IM,Hsobel,'replicate');
IMs2=imfilter(IM,Hsobel','replicate');

IMst=IMs+IMs2;
IMpt=IMp+IMp2;

IMstn=uint8(255.*ones(size(IM)))-IMst;
IMptn=uint8(255.*ones(size(IM)))-IMpt;

