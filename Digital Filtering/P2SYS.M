% 1 conjg poles pair
j=sqrt(-1);
fs=1024;
T=1./fs;

fp1=150;
sigp1=-10;

sp1=sigp1+j.*2.*pi.*fp1;
sp1c=sigp1-j.*2.*pi.*fp1;

zp1=exp(sp1.*T);
zp1c=exp(sp1c.*T);

% impulse responses

As = 1 ./ (sp1-sp1c);
%Az = 1 ./ (zp1-zp1c); same as s with scaling!

n=1:300;
ta=n.*T;
hs=As*exp(sp1.*n.*T)-As*exp(sp1c.*n.*T);
hz=As*(zp1.^(n-1))-As*(zp1c.^(n-1));

figure(1);
plot(ta,real(hz),'w.');
hold
plot(ta,real(hs),'y');
hold;

% system frequency responses

ifs=0:fs/1000:fs/2;
s=0+j.*2.*pi.*ifs;
z=exp(s.*T);

hw=1./((s-sp1).*(s-sp1c));
hwzc=( T.*(zp1-zp1c)./(sp1-sp1c) )./( (z-zp1).*(z-zp1c) );

figure(2);
plot(ifs,20.*log10(abs(hw)),'y');
hold;
plot(ifs,20.*log10(abs(hwzc)),'w.');
hold;

figure(3);
plot(ifs,57.296.*(angle(hw)),'y');
hold;
plot(ifs,57.296.*(angle(hwzc)),'w.');
hold;

% plot data to disk for Stanford Graphics

rowcol = [500 ; 5];
resp = [ifs; 20.*log10(abs(hw)); 57.296.*angle(hw); 20.*log10(abs(hwzc)); 57.296.*angle(hwzc)];
fd = fopen('z2resp.dat','wt');
fprintf(fd,'%d\t%d\n',rowcol);
fprintf(fd,'%e\t%e\t%e\t%e\t%e\n',resp);
fclose(fd);
