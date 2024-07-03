% treble shelf filter

fs=44100;
T=1/fs;

f0=5000;
G=.2;

a=(G-1)/(G+1);

j=sqrt(-1);
pi=acos(-1);
w0=2*pi*f0;

% treble filter
% single pole and zero
sz=-w0*(1-a);    % s-plane
sp=-w0*(1+a);
zz=exp(sz*T);   % mapped to z-plane
zp=exp(sp*T);

% evaluate frequency response
nfreqs=80;
fnyquist=fs/2;
F=[1:nfreqs].*(fnyquist/nfreqs);
sv=j.*2.*pi.*F;             % s-plane frequency
zv=exp(j.*2.*pi.*F./fs);        % z-plane frequency

scl=(1-zp)/(1-zz);          % make it unity gain at 0 Hz
Hs=G.*(sv - ones(size(sv)).*sz)./(sv - ones(size(sv)).*sp);

Hz=scl.*(zv - ones(size(zv)).*zz)./(zv - ones(size(zv)).*zp);

figure(1)
semilogx(F,abs(Hz));

% bass filter
% single pole and zero
sz=-w0*(1+a);    % s-plane
sp=-w0*(1-a);
zz=exp(sz*T);   % mapped to z-plane
zp=exp(sp*T);

% evaluate frequency response
nfreqs=80;
fnyquist=fs/2;
F=[1:nfreqs].*(fnyquist/nfreqs);
sv=j.*2.*pi.*F;             % s-plane frequency
zv=exp(j.*2.*pi.*F./fs);        % z-plane frequency

scl=(-1-zp)/(-1-zz);          % make it unity gain at fs/2 Hz
Hs=(sv - ones(size(sv)).*sz)./(sv - ones(size(sv)).*sp);

Hz=scl.*(zv - ones(size(zv)).*zz)./(zv - ones(size(zv)).*zp);

figure(2)
semilogx(F,abs(Hz));
