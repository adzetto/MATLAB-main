%upsampling of photos% open grayscale jpeg image

IM = imread('House Dec 2004 gs.jpg');

info = imfinfo('House Dec 2004 gs.jpg');

% NxN pixilation
[nr, nc]=size(IM);
% lets zoom in to see some detail
%nx1=round(.2*nr);
%nx2=round(.5*nr);
%ny1=round(.2*nc);
%ny2=round(.5*nc);

nx1=1;
nx2=nr;
ny1=1;
ny2=nc;

figure(1);
imshow(IM(nx1:nx2,ny1:ny2));
colormap(gray);
title('Original');

N=8;
nrdeci=fix(nr/N);
ncdeci=fix(nc/N);

nx1p=fix(nx1/N);    % for pixelated zoom in
if nx1p < 1
    nx1p=1;
end;
nx2p=fix(nx2/N);
ny1p=fix(ny1/N);
if ny1p < 1
    ny1p=1;
end;
ny2p=fix(ny2/N);

Hpix=(1/(N*N)).*ones([N,N]);
IMpix=zeros([nrdeci,ncdeci]);
scl=1/N^2;

IMzp=zeros(size(IM));   % this image will be oversampled N:1 of IMpix
for i=1:nrdeci
    for j=1:ncdeci
        sum=0;
        for ii=1:N
            for jj=1:N
                sum=sum+single(IM((i-1)*N+ii,(j-1)*N+jj));
            end;
        end;
        IMpix(i,j)=uint8(scl*sum); % do the math with precision, then 8-bit
        IMzp(N*i,N*j)=IMpix(i,j);
    end;
end;
    
% upsampling filters
F0=ones([8, 8]);            % 0th order
F1=conv2(F0,F0,'full');     % 1st order
F3=conv2(F1,F1,'full');     % 3rd order
F7=conv2(F3,F3,'full');     % 7th order

F0=F0./max(max(F0));
F1=F1./max(max(F1));
F3=F3./max(max(F3));
F7=F7./max(max(F7));

figure(2);
subplot(2,2,1),surf(F0);
colormap(white);
title('0th order');
subplot(2,2,2),surf(F1);
colormap(white);
title('1st order');
subplot(2,2,3),surf(F3);
colormap(white);
title('3rd order');
subplot(2,2,4),surf(F7);
colormap(white);
title('7th order');


figure(3);
IM0=filter2(F0,IMzp);
subplot(2,2,1),imagesc(IM0(nx1:nx2,ny1:ny2));
axis off;
colormap(gray);
title('0th order');

IM1=filter2(F1,IMzp);
subplot(2,2,2),imagesc(IM1(nx1:nx2,ny1:ny2));
axis off;
colormap(gray);
title('1st order');

IM3=filter2(F3,IMzp);
subplot(2,2,3),imagesc(IM3(nx1:nx2,ny1:ny2));
axis off;
colormap(gray);
title('3rd order');

IM7=filter2(F7,IMzp);
subplot(2,2,4),imagesc(IM7(nx1:nx2,ny1:ny2));
axis off;
colormap(gray);
title('7th order');

% OK try sharpening the updampled image
Hunsharp=fspecial('unsharp');
IMsharp=filter2(Hunsharp,IM3);
figure(4);
%image(IMupsharp(nx1:nx2,ny1:ny2));
imagesc(IMsharp);
axis off;
colormap(gray);
title('3rd order unsharp');

Hlap=fspecial('laplacian');
IMlap=filter2(Hlap,IM7);
IMlap=IM7-80.*IMlap;
figure(5);
%image(IMupsharp(nx1:nx2,ny1:ny2));
image(.07.*IMlap);
axis off;
colormap(gray);
title('7th order minus 80*LaPlacian');





