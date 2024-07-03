% echo and reverb
eval('clear all');

[y,fs,bits]=wavread('testing123.wav');
[npts nchan]=size(y);
%sound(y,fs);
%pause;

ysb=y;
t=0.075;
d=round(t*fs);
g=.8;

for n=d:npts,
    ysb(n,:)=y(n,:)+g.*y(n-d+1,:);
end;
%sound(ysb,fs);
wavwrite(ysb,fs,16,'slapback.wav');

yle=y;
t=0.500;
d=round(t*fs);
g=.9;

for n=d:npts,
    yle(n,:)=y(n,:)+g.*y(n-d+1,:);
end;
%sound(ysb,fs);
wavwrite(yle,fs,16,'longecho.wav');


spd1=.045; %spring delay
d1=round(spd1*fs);
fb1=0.7;

ysp1=zeros(size(y));
for n=d1:npts,
    ysp1(n,:)=y(n,:)+fb1.*ysp1(n-d1+1,:);
end;
%sound(ysp1,fs);
wavwrite(ysp1,fs,16,'springrev.wav');

spd2=.037; %spring delay
d2=round(spd2*fs);
fb2=0.65;

ysp2=zeros(size(y));
for n=d2:npts,
    ysp2(n,:)=y(n,:)+fb2.*ysp2(n-d2+1,:);
end;
y2sp=ysp1+ysp2;
%sound(y2sp,fs);
wavwrite(y2sp,fs,16,'platerev.wav');

spd3=.042; %spring delay
d3=round(spd3*fs);
fb3=0.6;

ysp3=y2sp;
for n=d3:npts,
    ysp3(n,1)=y2sp(n,1)+fb3.*ysp3(n-d3+1,1);    % process left only
    ysp3(n,2)=y2sp(n,2)+fb3.*ysp3(n-d3+10,2);    % process right only
end;
ysp3=ysp3.*.9; % clip control
sound(ysp3,fs);
wavwrite(ysp3,fs,16,'roomrev.wav');
