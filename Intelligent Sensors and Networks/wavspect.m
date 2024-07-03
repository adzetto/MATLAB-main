% make spectogram of bird wongs (wave files)

[file path]=uigetfile('*.wav','Input Wave File?');
pathfile=sprintf('%s%s',path,file);

[y fs bits]=wavread(pathfile);

[npts temp]=size(y);

nfft=1024;
%[YS F T]=specgram(y,nfft,fs,nfft,0);

figure(1);
%imagesc(T,F,20.*log10(abs(flipud(YS))));

[yf,f,t] = spectrogram(y,nfft,nfft/2,nfft,fs); 

imagesc(f,t,20.*log10(abs(yf')));
%surf(t,f,10*log10(abs(p)),'EdgeColor','none');   
%axis xy; 
%axis tight; 
colormap(bone); 
view(0,90);
ylabel('Time');
xlabel('Hz');

title('Spectrogram of Bird Song');

sound(y,fs);