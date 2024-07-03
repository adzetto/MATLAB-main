% parametric EQ in digital domain

fs=44100;

f0=9600;
Q=2;

G=[.1 .5 1 2 10].';
for k=1:5,
    if G(k) >= 1
        maz(k)=1;   % set az=1 if G>1 (peak)
        map(k)=maz(k)/G(k);
    else
        map(k)=1;   % set ap=1 if G<1 (dip)
        maz(k)=G(k)/map(k);
    end;
end;
az=maz.';
ap=map.';

% s-plane roots
j=sqrt(-1);
pi=acos(-1);
w0=2*pi*f0;

% numerator
sigz=-.5.*az./Q;
s1z=sigz.*w0+j.*w0.*sqrt(ones(size(G))-sigz.*sigz);
s2z=conj(s1z);

% denominator
sigp=-.5.*ap./Q;
s1p=sigp.*w0+j.*w0.*sqrt(ones(size(G))-sigp.*sigp);
s2p=conj(s1p);

% map to z-plane
T=1/fs;
z1z=exp(s1z.*T);
z2z=exp(s2z.*T);
z1p=exp(s1p.*T);
z2p=exp(s2p.*T);

% z-domain filter coefficients
bz=-z1z-z2z;
cz=z1z.*z2z;
bp=-z1p-z2p;
cp=z1p.*z2p;

% evaluate digital filter over nfreqs on unit circle
nfreqs=80;
fnyquist=fs/2;
F=[1:nfreqs].*(fnyquist/nfreqs);
zv=exp(j.*2.*pi.*F./fs);

% s-to-z plane scale factors
scl=( (z1p-z2p).*(s1z-s2z) )./( (s1p-s2p).*(z1z-z2z) );
scleye=zeros(size(scl*scl'));   % make identity matrix
for k=1:5,
    scleye(k,k)=scl(k);
end;

% generate frequency response

num=ones(size(G))*zv.^2 + bz*zv + cz*ones(size(zv));
den=ones(size(G))*zv.^2 + bp*zv + cp*ones(size(zv));
resp=scleye*num./den;

figure(1);
plot(F./1000,20.*log10(abs(resp(1,:))),'k-.');
ylabel('dB');
xlabel('kHz');
titbuf=sprintf('Parametric Responses Q=%.1f f_0=%.1f',Q,f0);
title(titbuf);
hold on;
plot(F./1000,20.*log10(abs(resp(2,:))),'k:');
plot(F./1000,20.*log10(abs(resp(3,:))),'k');
plot(F./1000,20.*log10(abs(resp(4,:))),'k-^');
plot(F./1000,20.*log10(abs(resp(5,:))),'k-o');
hold off;
legend('G=0.1 \alpha_z=0.1 \alpha_p=1.0',...
    'G=0.5 \alpha_z=0.5 \alpha_p=1.0',...
    'G=1.0 \alpha_z=1.0 \alpha_p=1.0',...
    'G=2.0 \alpha_z=1.0 \alpha_p=0.5',...
    'G=10.0 \alpha_z=1.0 \alpha_p=0.1');
