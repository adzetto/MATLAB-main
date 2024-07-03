npts=64;
fbins=1024;
trect=ones(size([1:npts]));

thann=zeros(size([1:npts]));
thamm=zeros(size([1:npts]));
tparz=zeros(size([1:npts]));
twelc=zeros(size([1:npts]));
texpn=zeros(size([1:npts]));
ttayl=zeros(size([1:npts]));
tscale=zeros(size([1:npts]));

kexp=-2.6;
for n=1:npts,
   thann(n) = .5*(1-cos( (2*pi*(n-1)/(npts-1)) ));
   thamm(n) = .5*(1.08-.92*cos( (2*pi*(n-1)/(npts-1)) ));
   tparz(n) = 1-abs( ((n-1)-.5*(npts-1))/(.5*(npts-1)) );
   twelc(n) = 1-( ((n-1)-.5*(npts-1))/(.5*(npts-1)) )^2;
   texpn(n) = sqrt(pi/3)*exp(-.5*3.43*(10^kexp)*(n-1-npts/2)^2);
   tscale(n) = (n-1)/npts;	
end;

figure(1);
plot(tscale,trect,'k');
hold on;
plot(tscale,thann,'k-o');
plot(tscale,tparz,'k:');
plot(tscale,twelc,'k--');
plot(tscale,texpn,'k-<');
plot(tscale,thamm,'k-s');
hold off;
line([0 0], [0 1],'Color','k');
line([tscale(npts) tscale(npts)], [1 0],'Color','k');
xlabel('n/N');
axis([-.1 1.1 0 1.1]);
legend('Rect','Hann','Parz','Welc','Expn','Hamm');
ylabel('W[n]');

fscale=zeros(size([1:fbins]));
for n=1:fbins,
   fscale(n) = (n-.5*fbins-1)/(fbins/npts);
end

frect=fft(trect,fbins);
frect=(1/npts).*fftshift(frect);
frect=abs(frect);

fhann=fft(thann,fbins);
fhann=(1/npts).*fftshift(fhann);
nhann=1/abs(fhann(1+fbins/2));
fhann=nhann.*abs(fhann);

fparz=fft(tparz,fbins);
fparz=(1/npts).*fftshift(fparz);
nparz=1/abs(fparz(1+fbins/2));
fparz=nparz.*abs(fparz);

fwelc=fft(twelc,fbins);
fwelc=(1/npts).*fftshift(fwelc);
nwelc=1/abs(fwelc(1+fbins/2));
fwelc=nwelc.*abs(fwelc);

fexpn=fft(texpn,fbins);
fexpn=(1/npts).*fftshift(fexpn);
nexpn=1/abs(fexpn(1+fbins/2));
fexpn=nexpn.*abs(fexpn);

fhamm=fft(thamm,fbins);
fhamm=(1/npts).*fftshift(fhamm);
nhamm=1/abs(fhamm(1+fbins/2));
fhamm=nhamm.*abs(fhamm);

flow=.5*fbins-4*fbins/npts;
fhi=.5*fbins+4*fbins/npts;

figure(2);
plot(fscale(flow:fhi),frect(flow:fhi),'k');
hold on;
plot(fscale(flow:fhi),fhann(flow:fhi),'k-o');
plot(fscale(flow:fhi),fparz(flow:fhi),'k:');
plot(fscale(flow:fhi),fwelc(flow:fhi),'k--');
plot(fscale(flow:fhi),fexpn(flow:fhi),'k-<');
plot(fscale(flow:fhi),fhamm(flow:fhi),'k-s');
hold off;
xlabel('Relative Bin');
ylabel('|Y(f)|');
axis([-4 4 0 1.1]);
%legend('Rectangular','Hanning');
legend('Rect','Hann','Parz','Welc','Expn','Hamm');

rectnb=1/(sum(trect)/npts)
rectbb=1/sqrt(sum(trect.^2)/npts)

welcnb=1/(sum(twelc)/npts)
welcbb=1/sqrt(sum(twelc.^2)/npts)

parznb=1/(sum(tparz)/npts)
parzbb=1/sqrt(sum(tparz.^2)/npts)

hammnb=1/(sum(thamm)/npts)
hammbb=1/sqrt(sum(thamm.^2)/npts)

hannnb=1/(sum(thann)/npts)
hannbb=1/sqrt(sum(thann.^2)/npts)

expnnb=1/(sum(texpn)/npts)
expnbb=1/sqrt(sum(texpn.^2)/npts)

