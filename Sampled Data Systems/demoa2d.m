% Matlab m-file for Figures 1.2 and 1.3 A2D-Demo

fs=1000;		% sample rate
Ts=1/fs;		% sample time interval

fs_analog=10000;	% "our" display sample rate (analog signal points)
npts_analog=200;	% number of analog display points
T_analog=1/fs_analog;	% "our" display sample interval

f0=950;             % use 75 Hz for Fig 1.3 and 950 Hz for Fig 1.4
Tstop=0.015;		% show 15 msecs of data

Ta=0:T_analog:Tstop;	% analog "samples"
Td=0:Ts:Tstop;		% digital samples

ya=zeros(size(Ta));	% zero out data vectors same length as time
yd=zeros(size(Td));

w0=2*pi*f0;

ya=cos(w0.*Ta);		% note scalar by vector multiply (.*) gives vector in 
                    % the cosine argument and a vector in the output ya

yd=cos(w0.*Td);

figure(1);          % initialize a new figure window for plotting

plot(Ta,ya,'k');    % plot in black
hold on;            % keep the current plot and add another layer
plot(Td,yd,'k*');   % plot in black "*"
hold off;           % return figure to normal state
xlabel('Seconds');

