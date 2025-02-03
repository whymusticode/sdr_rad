
clear

c = 3e8; % speed of light m/s

fc = .7e9; % center frequency
fs = 6e6; % sample frequency 

r = 10e3; % target range m
v = -10; % target velocity at start of CPI
a = 0; % target acceleration m/s^2
SNR = 1e-0; % target signal relative to noise floor in voltage amplitude

rMax = 40e3; % maximum range
vMax = 20; % maximum velocity m/s

cpi = 6;

n = round(fs*cpi); % number of samples in integration time
df = fc*v/c; % doppler HZ
df2 = fc*(v+cpi*a)/c; % doppler HZ

% create transmit signal
sig = randn(n,1) + 1i*randn(n,1);

noise = randn(n,1) + 1i*randn(n,1);
t = (1:n)'/n * cpi;
dShift = exp(1i*2*pi*  (df*t + t.^2/2*(df2-df)/cpi)  );
rShift =[zeros(round(r/c*fs),1);sig(1:end-round(r/c*fs))];
RX = rShift.*dShift*SNR + noise + 0*sig; % recieved signal
sig = single(sig);
RX = single(RX); % single precision is faster
%% Signal Processing

p.rMax = rMax; p.fc = fc; p.fs = fs; p.vMax = vMax;
% profile on
tic
[RANGES,VELS,z2] = sledgeBatch(sig,RX,p);
% toc
% profile viewer
t1 = toc;
% return

% addpath(genpath(['C:\Users\mi25980\Desktop\EA\0LD\',...
%   'analyzeWaveforms\Sledgehamer']))
% radar.Ptx = 1;
% p.l = c/p.fc;
% radar.Dmax = max(VELS)/2;
% radar.Rmax = p.rMax;
% IQRX.full = RX;
% IQTX = sig;
% settings.filter_clutter = false;
% settings.range_taper = false;
% settings.doppler_taper = true;
% settings.STAR = true;
% pulse.f0  = p.fc - p.fs/2;
% pulse.bw = p.fs;
% tic
% [d_RDI, IQfilt, RDI]= sledgehammer_processing(radar, p, pulse, ...
%   t, IQTX, true(size(t)), false(size(t)),12, IQRX, settings);
% t2 = toc;
% slR = abs(RDI.full(2:size(z2,2)+1,:)');
% snr1 = 20*log10(max(abs(z2(:)))/median(abs(z2(:))))
% snr2 = max(abs(slR(:)))/median(abs(slR(:)));
% disp(['speed improvement factor of ',num2str(round(t2/t1))])
% disp(['snr is [new current]:  ', num2str(round([snr1,snr2]))])



figure(2)
az2 = abs(z2);
imagesc(RANGES/1e3,VELS,20*log10(az2/median(az2(:))))
xlabel('Range Difference of Arrival (km)')
ylabel('Doppler shift (m/s)')
title('fast algorithm')
colorbar
caxis([0,30])
% xlim([40,60])
% ylim([180,240])
colormap(jet)
font
% subplot(1,2,2)
% imagesc(RANGES/1e3,-d_RDI*2,20*log10(slR))
% % xlim([40,60])
% % ylim([180,240])
% title('current sledgehammer')
% colorbar
% caxis([5,50])
% colormap(jet)


