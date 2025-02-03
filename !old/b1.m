
% this stores rx data (radio address '01')
% (master oscillator card, with resistor)

clear
load params
sdrRx.RadioID = 'usb:1';
saveName = './data/tx';
saveDataLoop(sdrRx,saveName,0,gLim)

% 
% clear
% 
% time = 100; % seconds of data to collect
% 
% fc = 96.6e6; 
% fs = 2.4e6;
% N = 375e3;
% g = 20; % amplification at the beginning (db?)
% 
% sdrRx = comm.SDRRTLReceiver('CenterFrequency', fc, 'SampleRate',fs, ...
%   'SamplesPerFrame',N,'OutputDataType',  'single'...
%   ,'TunerGain',g,'EnableTunerAGC',  false);%,'TunerGain',g,'EnableTunerAGC',  false
% sdrTx = comm.SDRRTLReceiver('CenterFrequency', fc, 'SampleRate',fs,...
%   'SamplesPerFrame',N,'RadioAddress','01','OutputDataType', ...
%   'single','TunerGain',g,'EnableTunerAGC',  false);%,'TunerGain',g,'EnableTunerAGC',  false
% CLOCK = clock;
% 
% 
% % LAG = 0;
% % setLag = 0;
% % thresh = 5000;
% 
% steps = round(time*fs/N);
% rx = zeros(steps*N,1);
% tx = rx;
% lost = zeros(steps,2);
% for i = 1:steps
%   [rx((1:N)+N*(i-1)), sdrRx,lost(i,1)] = myStep(sdrRx);
%   [tx((1:N)+N*(i-1)), sdrTx,lost(i,2)] = myStep(sdrTx);
% end
% 
% [t,r] = alignFun(rx,tx,10*fs);



% sum(lost(:))
% lag = 0*fs;
% buf = 1e8;
% integ = 20*fs;
% RX = fft(rx((1:integ)+buf));
% TX = fft(conj(flipud(tx((1:integ)+buf))));
% res = ifft(RX.*TX);
% short = smpl(abs(res),1e3);

% [~,coarseLag] = max(abs(res));

% plot(short)












%   i2 = floor(i/2);
%   if ~setLag
%     figure(506)
%     subplot(2,1,1)
%     [lag,~,~,MAX] = getLag(rx((1:N)+N*i2),tx(end+1-N:end),1);
%     if MAX > thresh
%       LAG = lag + (i-1-i2)*N;
%       setLag = 1;
%     end
%     subplot(2,1,2)
%     [lag,~,~,MAX] = getLag(rx(end+1-N:end),tx((1:N)+N*i2),1);
%     if MAX >thresh
%       LAG = lag - (i-1-i2)*N;
%       setLag = 1;
%     end
%     continue
%   end
%   
%   if LAG > 0;  D1 = rx(end-N+1:end); D2 = tx(end-N+1-LAG:end-LAG);
%   else;        D2 = tx(end-N+1:end); D1 = rx(end-N+1+LAG:end+LAG);
%   end
%   
%   figure(506)
%   subplot(1,1,1)
%   
%   [lag,~,~,MAX] = getLag(D1,D2,1);
  




% 
% getLag([1,0,0,0]',[0,0,1,0]',0)
% 
% 
% plot(real(tx))

%
% % tic
% getLag(data(:,1),data(:,2),1)
% % toc
%
% length(unique((real(tmp1))))
%
% [vals,lags] = xcorr(data(:,1),data(:,2));
% plot(lags,abs(vals))


% f = randn(100,1);
% [vals,lags] = xcorr(f,f);
% plot(lags,abs(vals))