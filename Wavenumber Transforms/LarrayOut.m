% line array beam output for 3 sources

thetad = 1:360;
theta = pi.*thetad./180;

M=16;
M=32;
doverl=.4;

S=zeros([4,1]);
ts=zeros([4,1]);

S(1)=10;
ts(1)=45;		% source 1 @ 65 deg level 30
S(2)=20;
ts(2)=85;		% source 2 @ 85 deg level 20
S(3)=25;
ts(3)=105;	% source 3 @ 105 deg level 25
S(4)=15;
ts(4)=145;	% source 4 @ 145 deg level 15

% here we sum up sources at the given directions of arrival (DOA)
y=zeros([1,M]);	% array data vector
% 4 sources, 16 array elements to sum source for
for n=1:4,
  steer = +j*2*pi*doverl*cos(pi*ts(n)/180); % phase offset for source
  for m=1:M,
    y(m) = y(m) + S(n)*exp(steer*(m-1));    % phase for mth element
  end;
end;

Ah=hanning(M);			% scaled Hanning shading
Ar=ones([M,1]);       % in case you want to see no Hanning window
Aw=zeros([M,1]);
Awr=zeros([M,1]);

scale = 0;
for m=1:M,
 scale = scale + Ah(m);
end;
Ah = Ah.*M./scale;		% amplitude scaled

resp = zeros(size(theta));
respr = zeros(size(theta));

% now we make a steering vector for each degree and "scan" 1 to 180
for n=1:180,
  steer = -j*2*pi*doverl*cos(pi*n/180);
  for m=1:M,
    Aw(m) = Ah(m)*exp(steer*(m-1)); % make steering vector
    Awr(m) = Ar(m)*exp(steer*(m-1)); % make steering vector (rect window)
  end;			
  resp(n)=y*Aw./M; % this is the array output for the nth look direction
  respr(n)=y*Awr./M; % this is the array output - rectangular window
end;

figure(1);
plot(thetad(1:180),abs(resp(1:180)),'k',thetad(1:180),abs(respr(1:180)),'k:');
ylabel('|s| / M');
xlabel('scanned theta');

titbuf=sprintf('M=%d-Element Windowed Line Array Output',M);
title(titbuf);
text(30,20,'d/\lambda = 0.4');
legend('Hanning','Rectangular');
axis([0 180 0 30]);
