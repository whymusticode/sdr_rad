function [r,v,rdi] = sledgeOld(rx,tx,p)
% clear


% p.fs = .2e6; % sample frequency
% p.fc = 1e8; % center frequency
% p.rMax = 1e5; % maximum range
% p.vMax = 600; % maximum doppler


c = 3e8;

%% Create real data
n = length(rx); % number of samples 
cpi = n/p.fs;
dMax = p.vMax/c*p.fc;
rBins = round(p.rMax/c*p.fs); % number of bins to go
DOPS = -dMax:1/(4*cpi):dMax;
dBins =  length(DOPS);

rx = gpuArray(rx);
tx = gpuArray(tx);

% DSHIFT = gpuArray(exp(1i*(1:n)')*gpuArray(DOPS*2*pi*cpi/n)); 
DSHIFT = gpuArray(exp(1i*(1:n)'*DOPS*2*pi*cpi/n)); 

% large matrix gpuArray(...
rep = bsxfun(@times,DSHIFT,tx);
x = repmat(rx,[1,dBins]);

%% for windowing
% WIN = repmat(hamming(n),[1,dBins]);
% rep = rep.*WIN;

rdi = ifft(fft(x,2*n).*fft(conj(flipud(rep)), 2*n));


% rdi = rdi(n+1:n+rBins,:)';
% r = (1:rBins)*c/2/p.fs;
buf = 0;
rdi = rdi(n-buf:n+rBins,:)';
r = (-buf:rBins)*c/p.fs;
v = DOPS*c/p.fc;




% figure(1)
% aRDI = abs(rdi);
% hold off
% imagesc(r/1e3,v,20*log10(aRDI/median(aRDI(:))))
% xlabel('Range Difference of Arrival (km)')
% ylabel('Doppler shift (m/s)')
% colorbar
% % caxis([0,20])
% colormap(jet)


%% This shows how a time shift is created by a phase shift in freq
%  t=0:0.001:0.1-0.001;
%  Fs = 1e3;
%  freq1 = 100;
%  x1=sinc(2*pi*freq1*t);
%  Delay=1.2;
%  yp = fft(x1);
%  yp = yp(1:length(x1)/2+1);
%  f = 0:Fs/length(x1):500;
%  yp = yp.*exp(-1i*2*pi*f*Delay*(1/Fs));
%  yp = [yp conj(fliplr(yp(2:end-1)))];
%  rep = ifft(yp,'symmetric');
%  hold off
%  plot(t(1:100),sinc(2*pi*freq1*(t-Delay/Fs)),'b:');
%  hold on;
%  plot(t(1:100),rep(1:100),'r');

