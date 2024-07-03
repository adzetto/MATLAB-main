% LS fit for A, B, and C frequency weighting curves

fd = [.01995 .05012 .1 .1995 1 1.995 3.981 6.31 10 20];		% freqs in kHz
logfd = log10(fd);

As = [-50.5 -30.2 -19.1 -10.9 0 1.2 1 -.1 -2.5 -9.3]';	% spl's below 55 dB
Bs = [-24.2 -11.6 -5.6 -2 0 -.1 -.7 -1.9 -4.3 -11.1]';	% spl's between 55 dB and 85 dB
Cs = [-6.2 -1.3 -.3 0 0 -.2 -.8 -2 -4.4 -11.2]';			% spl's greater than 85 dB ref: Beranek&Ver, p15

M = 5;				% model order
N = 10;				% num of frequencies
Fbar = zeros([N M]);

for m=1:M,
  for n=1:N
    Fbar(n,m) = logfd(n)^(m-1);	% Vandermode matrix of log freqs
  end;
end;

FFIF = (Fbar'*Fbar)^(-1)*Fbar';	% least-squares solution

Ha = FFIF*As;
Hb = FFIF*Bs;
Hc = FFIF*Cs;

nfpts = 1000;
df = fd(N)/nfpts;
for n=1:nfpts,
  fa(n) = df*n;					% frequency axis (linear scale)
  logfa(n) = log10(fa(n));
  sa(n)=0;
  sb(n)=0;
  sc(n)=0;						% sa, ab, sc are the frequency weightings
  for m=1:M,
    sa(n)=sa(n)+Ha(m)*logfa(n)^(m-1);
    sb(n)=sb(n)+Hb(m)*logfa(n)^(m-1);
    sc(n)=sc(n)+Hc(m)*logfa(n)^(m-1);
  end;
end;
								% make semilog plots for A, B, and C wtgs
figure(1);
semilogx(fa,sa,'k');
hold on;
semilogx(fa,sb,'k');
semilogx(fa,sc,'k');
semilogx(fd, As,'k*');
semilogx(fd, Bs,'ko');
semilogx(fd, Cs,'k+');
hold off;
								% label axis, add legend & title
xlabel('kHz');
ylabel('dB Weighting');
title('dB Adjustments for A, B, and C Weightings');
legend('A[*] SPL < 55dB', 'B[o] 55dB < SPL < 85dB', 'C[+] SPL > 85dB');

								% write weights out to disk

FP=fopen('abcwts.dat','wt');	% ASCII text file format
fprintf(FP,'%d\n',M);			% model order
for m=1:M,
  fprintf(FP,'%+e\t%+e\t%+e\n',Ha(m),Hb(m),Hc(m));
end;
fclose(FP);
