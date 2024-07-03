nfpts=256;
nbuff=256;
npts=nbuff*nfpts;
fs=128;

% 1st sine wave
a1=12;
f1=10;
th1=45*pi/180;

% 2nd sine wave
a2=10;
f2=14;
th2=170*pi/180;

% nonlinearity and noise levels
eps=.035;
nz=0.01;
alpha=10;
beta=0.05;
%titbuf=sprintf('Nz=%.2f \\epsilon=%.3f',nz,eps);
titbuf=sprintf('Nz=%.2f \\alpha=%.1f \\beta=%.3f',nz,alpha,beta);

w1=2*pi*f1/fs;
w2=2*pi*f2/fs;

% generate the clean, distorted, and a fake input with the same freqs as
% the distorted
for n=1:npts,
   x(n)=a1*cos(w1*n+th1)+a2*cos(w2*n+th2)+nz*randn;
   z(n)=x(n)+eps*x(n)^2;
   %z(n)=alpha*atan(beta*x(n));
   zfake(n)=a1*cos(w1*n+th1)+a1*cos(2*w1*n+2*th1)+.5*(a1+a2)*cos((w1-w2)*n+(th1-th2))+ ... 
       .5*(a1+a2)*cos((w1+w2)*n+(th1+th2))+ a2*cos(w2*n+th2)+ ...
        a2*cos(2*w2*n+2*th2) + nz*randn;
end
win=hanning(nfpts);
for n=1:nbuff,
    fbuff=fft(win'.*x((n-1)*nfpts+1:n*nfpts),nfpts);
    fx((n-1)*nfpts+1:n*nfpts)=fbuff;
    fbuff=fft(win'.*z((n-1)*nfpts+1:n*nfpts),nfpts);
    fz((n-1)*nfpts+1:n*nfpts)=fbuff;
    fbuff=fft(win'.*zfake((n-1)*nfpts+1:n*nfpts),nfpts);
    fzfake((n-1)*nfpts+1:n*nfpts)=fbuff;
end;

% lets compare time waveforms
np=100;
T=(1:np)./fs;

figure(1);
%plot(T,x(1:np),'k',T,z(1:np),'k:',T,zfake(1:np),'k-.');
%legend('Original','Distorted','Linear Sines');
plot(T,x(1:np),'k',T,z(1:np),'k:');
legend('Original','Distorted');
xlabel('Sec');
title(titbuf);

% Let's compare Power spectra
Px=zeros(1,nfpts/2);
Pz=zeros(1,nfpts/2);
Pzfake=zeros(1,nfpts/2);

for k=1:nbuff,
    Px=Px+abs(fx(1:nfpts/2)).^2;
    Pz=Pz+abs(fz(1:nfpts/2)).^2;
    Pzfake=Pzfake+abs(fzfake(1:nfpts/2)).^2;
end;
Px=Px./nbuff;
Pz=Pz./nbuff;
Pzfake=Pzfake./nbuff;
% mult fft by (2/nfpts)^2 for NPSD
Px=4.*Px./(nfpts^2);
Pz=4.*Pz./(nfpts^2);
Pzfake=4.*Px./(nfpts^2);
F=(1:nfpts/2).*fs./nfpts;

figure(2);
%plot(F,10.*log10(Px),'k',F,10.*log10(Pz),'k:',F,10.*log10(Pzfake),'k-.');
%legend('Original','Distorted','Linear Sines');
plot(F,10.*log10(Px),'k',F,10.*log10(Pz),'k:');
legend('Original','Distorted');
xlabel('Hz');
ylabel('NPSD dB');
title(titbuf);

% now lets compute the bispectrum

bispx=zeros(nfpts/2,nfpts/2);
bispz=zeros(nfpts/2,nfpts/2);
bispzfake=zeros(nfpts/2,nfpts/2);

for n=1:nbuff,
    for k=1:nfpts/4,
        for kk=1:nfpts/4,
            bispx(k,kk)=bispx(k,kk)+fx((n-1)*nfpts+k)*fx((n-1)*nfpts+kk)*...
                conj(fx((n-1)*nfpts+k+kk));
            bispz(k,kk)=bispz(k,kk)+fz((n-1)*nfpts+k)*fz((n-1)*nfpts+kk)*...
                conj(fz((n-1)*nfpts+k+kk));
            bispzfake(k,kk)=bispzfake(k,kk)+fzfake((n-1)*nfpts+k)*fzfake((n-1)*nfpts+kk)*...
                conj(fzfake((n-1)*nfpts+k+kk));
        end;
    end;
end;

    for k=1:nfpts/4,
        for kk=1:nfpts/4,
            bispcox(k,kk)=bispx(k,kk)/sqrt(Px(k)*Px(kk)*Px(k+kk));
            bispcoz(k,kk)=bispz(k,kk)/sqrt(Pz(k)*Pz(kk)*Pz(k+kk));
            bispcozfake(k,kk)=bispzfake(k,kk)/sqrt(Pzfake(k)*Pzfake(kk)*Pzfake(k+kk));
        end;
    end;

[Fx,Fy]=meshgrid(F(1:nfpts/4),F(1:nfpts/4));
figure(3);
colormap(jet);
%pcolor(Fx,Fy,10.*log10(abs(bispx(1:nfpts/4,1:nfpts/4))));
pcolor(Fx,Fy,unwrap(angle(bispx(1:nfpts/4,1:nfpts/4))));
shading interp;
caxis([0,100]);
colorbar;

figure(4);
colormap(jet);
%pcolor(Fx,Fy,10.*log10(abs(bispz(1:nfpts/4,1:nfpts/4))));
pcolor(Fx,Fy,unwrap(angle(bispz(1:nfpts/4,1:nfpts/4))));
shading interp;
caxis([0,100]);
colorbar;

figure(5);
colormap(jet);
%pcolor(Fx,Fy,10.*log10(abs(bispzfake(1:nfpts/4,1:nfpts/4))));
pcolor(Fx,Fy,unwrap(angle(bispzfake(1:nfpts/4,1:nfpts/4))));
shading interp;
caxis([0,100]);
colorbar;

% let's look at 1 axis
figure(6);
%plot(F(1:nfpts/4),10.*log10(abs(bispx(1,1:nfpts/4))),'k',...
%    F(1:nfpts/4),10.*log10(abs(bispz(1,1:nfpts/4))),'k:',...
%    F(1:nfpts/4),10.*log10(abs(bispzfake(1,1:nfpts/4))),'k-.');
%legend('Original','Distorted','Linear Sines');
plot(F(1:nfpts/4),10.*log10(abs(bispx(1,1:nfpts/4))),'k',...
    F(1:nfpts/4),10.*log10(abs(bispz(1,1:nfpts/4))),'k:');
legend('Original','Distorted');
xlabel('Hz');
ylabel('10log_10|C_3^x(0,f_2)|');
title(titbuf);

% let's look at 1 axis of bispec coher
figure(7);
%plot(F(1:nfpts/4),10.*log10(abs(bispx(1,1:nfpts/4))),'k',...
%    F(1:nfpts/4),10.*log10(abs(bispz(1,1:nfpts/4))),'k:',...
%    F(1:nfpts/4),10.*log10(abs(bispzfake(1,1:nfpts/4))),'k-.');
%legend('Original','Distorted','Linear Sines');
plot(F(1:nfpts/4),10.*log10(abs(bispcox(1,1:nfpts/4))),'k',...
    F(1:nfpts/4),10.*log10(abs(bispcoz(1,1:nfpts/4))),'k:');
legend('Original','Distorted');
xlabel('Hz');
ylabel('10log_10|C_3^x(0,f_2)|');
title(titbuf);
