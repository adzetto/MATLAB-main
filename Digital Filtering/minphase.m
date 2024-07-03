% 1 conjg poles pair
% 1 conjg zeros pair

j=sqrt(-1);
fs=1000;
T=1./fs;

% original system
fz1=300;
sigz1=+60;
sz1=sigz1+j.*2.*pi.*fz1;
sz1c=sigz1-j.*2.*pi.*fz1;
zz1=exp(sz1.*T);
zz1c=exp(sz1c.*T);

% all pass system
zp1=zz1;
zp1c=zz1c;

% min phase system
zzm1=1./conj(zz1);
zzm1c=1./conj(zz1c);
apscale=1/abs(zzm1)^2;
%apscale=1;

% system frequency responses

ifs=0:fs/100:fs/2;
s=0+j.*2.*pi.*ifs;
z=exp(s.*T);

hw = (1.-zz1./z).*(1-zz1c./z);
hap = apscale.*((z-zzm1).*(z-zzm1c))./((z-zp1).*(z-zp1c));
hm = hw.*hap;

figure(1);
subplot(2,1,2),plot(ifs,20.*log10(abs(hw)),'k',...
    ifs,20.*log10(abs(hap)),'k-.',ifs,20.*log10(abs(hm)),'k-o');
legend('original system','all pass system','min phase system');
xlabel('Hz');
ylabel('dB');
subplot(2,1,1),plot(ifs,57.296.*(unwrap(angle(hw))),'k',ifs,...
    57.296.*(unwrap(angle(hap))),'k-.',ifs,57.296.*(unwrap(angle(hm))),'k-o');
ylabel('Degrees');
xlabel('Hz');
legend('original system','all pass system','min phase system');
