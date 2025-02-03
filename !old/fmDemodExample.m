
clear
fs = 10000; % sample
fc = 200;  % carrier
fDev = 50; % deviation

t = (0:1/fs:0.2)';

x = sin(2*pi*30*t)+2*sin(2*pi*60*t);


y = fmmod(x,fc,fs,fDev);
yM = y + randn(size(y))*.1;
xEst = fmdemod(yM,fc,fs,fDev); % Demodulate both channels.
yEst = fmmod(xEst,fc,fs,fDev);

plot(t,yM,'.',t,yEst)
plot(t,x,'.',t,xEst)


[val,lag] = xcorr(yM-mean(yM),yEst-mean(yEst));
% [val,lag] = xcorr(x-mean(x),z-mean(z));
plot(lag,val)
font


% plot(t,x,t,z)
