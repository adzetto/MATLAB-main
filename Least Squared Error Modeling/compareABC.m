% test spectrum for A, B, and C weighting

FP=fopen('abcwts.dat','rt');	% open as ASCII text
M=fscanf(FP,'%d',1);
Ha=zeros([M 1]);
Hb=zeros([M 1]);
Hc=zeros([M 1]);				% mx1 column vectors
for m=1:M,
  Ha(m)=fscanf(FP,'%e',1);
  Hb(m)=fscanf(FP,'%e',1);
  Hc(m)=fscanf(FP,'%e',1);
end;
fclose(FP);

fd = [.01995 .05012 .1 .1995 1 1.995 3.981 6.31 10 20];		% freqs in kHz
As = [-50.5 -30.2 -19.1 -10.9 0 1.2 1 -.1 -2.5 -9.3]';	% spl's below 55 dB
Bs = [-24.2 -11.6 -5.6 -2 0 -.1 -.7 -1.9 -4.3 -11.1]';	% spl's between 55 dB and 85 dB
Cs = [-6.2 -1.3 -.3 0 0 -.2 -.8 -2 -4.4 -11.2]';			% spl's greater than 85 dB ref: Beranek&Ver, p15


% generate npts of test data - sinusoids in white noise

npts=1024;
y=zeros([1 npts]);
fs=20000;
T=1/fs;
f=fs/npts:fs/npts:fs;

% create A, B, and C, weighting filters

logf=log10(f./1000);
sa=zeros([1 npts]);
sb=zeros([1 npts]);
sc=zeros([1 npts]);

for n=1:npts,
  sa(n)=0;
  sb(n)=0;
  sc(n)=0;						% sa, ab, sc are the frequency weightings
  for m=1:M,					% in dB
    sa(n)=sa(n)+Ha(m)*logf(n)^(m-1);
    sb(n)=sb(n)+Hb(m)*logf(n)^(m-1);
    sc(n)=sc(n)+Hc(m)*logf(n)^(m-1);
  end;
end;

fd1k=fd*1000; % for overlay graphics
nfd1k=round(npts.*fd1k./fs);

figure(1);
semilogx(f,sa,'k');
hold on;
semilogx(f,sb,'k:');
semilogx(f,sc,'k-');
semilogx(fd1k, As,'k*');
semilogx(fd1k, Bs,'ko');
semilogx(fd1k, Cs,'k+');
hold off;
title('LS-fit ABC Curves');
ylabel('dB');
xlabel('Hz');
legend('A[*]', 'B[o]', 'C[+]');

% now let's compare to the ANSI S1.4-1983

wa=zeros([1 npts]);
wb=zeros([1 npts]);
wc=zeros([1 npts]);

k1=2.242881*10^16;
k2=1.025119;
k3=1.562339;
f1=20.598997;
f2=107.65265;
f3=737.86223;
f4=12194.22;
f5=158.48932;

for n=1:npts,
  wc(n) = 10.*log10( (k1*f(n)^4)/( ((f(n)^2+f1^2)^2)*((f(n)^2+f4^2)^2) ) );
  wb(n) = 10.*log10( (k2*f(n)^2)/(f(n)^2+f5^2) ) + wc(n);
  wa(n) = 10.*log10( (k3*f(n)^4)/((f(n)^2+f2^2)*(f(n)^2+f3^2)) ) + wc(n);
end;

figure(2);
semilogx(f,wa,'k');
hold on;
semilogx(f,wb,'k');
semilogx(f,wc,'k');
semilogx(fd1k, As,'k*');
semilogx(fd1k, Bs,'ko');
semilogx(fd1k, Cs,'k+');
hold off;
title('ANSI S1.1983 ABC Curves');
ylabel('dB');
xlabel('Hz');
legend('A[*]', 'B[o]', 'C[+]');

figure(3);
semilogx(f,(wa-sa),'k');
hold on;
semilogx(f,(wb-sb),'k');
semilogx(f,(wc-sc),'k');
semilogx(fd1k,wa(nfd1k)-sa(nfd1k) ,'k*');
semilogx(fd1k,wb(nfd1k)-sb(nfd1k) ,'ko');
semilogx(fd1k,wc(nfd1k)-sc(nfd1k) ,'k+');
hold off;
title('ANSI Minus LS-fit ABC Curves');
ylabel('dB');
xlabel('Hz');
legend('A[*]', 'B[o]', 'C[+]');
axis([10 100000 -.6 .6]);


