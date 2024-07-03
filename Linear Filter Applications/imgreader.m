% open grayscale jpeg image

IM = imread('House Dec 2004 gs.jpg');

info = imfinfo('House Dec 2004 gs.jpg');
figure(1);
imshow(IM);
colormap(gray);

% NxN pixilation
[nr, nc]=size(IM);
N=16;
nrdeci=fix(nr/N);
ncdeci=fix(nc/N);
Hpix=(1/(N*N)).*ones([N,N]);
IMpix=zeros([nrdeci,ncdeci]);
scl=1/N^2;
for i=1:nrdeci
    for j=1:ncdeci
        sum=0;
        for ii=1:N
            for jj=1:N
                sum=sum+single(IM((i-1)*N+ii,(j-1)*N+jj));
            end;
        end;
        IMpix(i,j)=uint8(scl*sum); % do the math with precision, then 8-bit
    end;
end;
    
figure(2);
imagesc(IMpix);
axis off;
colormap(gray);

% boxcar low pass filter
Hbc=ones([N,N]);
IMbc=imfilter(IM,Hbc);
figure(3);
imshow(IMbc)



Hunsharp=fspecial('unsharp');
IMs=imfilter(IM,Hunsharp,'replicate');

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

