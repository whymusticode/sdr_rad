


% 87.9 - 107.9


reset(gpuDevice)
clear

time = 120; % seconds of data to collect

p.fc = 97.3*1e6; %389+6*rfChannel
p.fs = 2.4e6;
N = 375e3;
g = 30; % amplification at the beginning (db?)
useGain = 1;


p.rMax = 1e5; % maximum range
p.vMax = 600; % maximum doppler
p.cpi = 2; % seconds to integrate


if useGain
sdrRx = comm.SDRRTLReceiver('CenterFrequency', p.fc, 'SampleRate',p.fs, ...
  'SamplesPerFrame',N,'OutputDataType',  'single'...
  ,'EnableTunerAGC',  false,'TunerGain',g);
sdrTx = comm.SDRRTLReceiver('CenterFrequency', p.fc, 'SampleRate',p.fs,...
  'SamplesPerFrame',N,'RadioAddress','01','OutputDataType', ...
  'single','EnableTunerAGC',  false,'TunerGain',g);
else
  sdrRx = comm.SDRRTLReceiver('CenterFrequency', p.fc, 'SampleRate',p.fs, ...
  'SamplesPerFrame',N,'OutputDataType',  'single'...
  );
sdrTx = comm.SDRRTLReceiver('CenterFrequency', p.fc, 'SampleRate',p.fs,...
  'SamplesPerFrame',N,'RadioAddress','01','OutputDataType', ...
  'single');
end

CLOCK = clock;

steps = round(time*p.fs/N);
rx = zeros(steps*N,1);
tx = rx;
lost = zeros(steps,2);
for i = 1:steps
  [rx((1:N)+N*(i-1)), sdrRx,lost(i,1),check1] = myStep(sdrRx);
  [tx((1:N)+N*(i-1)), sdrTx,lost(i,2),check2] = myStep(sdrTx);
end
disp(['lost this many samples: ',num2str(sum(lost(:)))])



tic
% rx = gather(rx);
% tx = gather(tx);
if check2 > check1 % use the sdr using more dynamic range as recieve 
  [txa,rxa] = alignFun(rx,tx);
else
  [rxa,txa] = alignFun(rx,tx);
end
toc

% alignFun(rxa(1:10:2e4),txa(1:10:2e4));
% rdi = LSsledge(txa(1:2e4),rxa(1:2e4),40,50);
% 
% 
% figure(461)
% imagesc(20*log10(abs(rdi')))
% colorbar


clear rx tx

rxa = int8(128*rxa);
txa = int8(128*txa);
rxa = rxa(1:1.5e8);
txa = txa(1:1.5e8);

save('dataCollect')

return 
% p.fs = 2.4e6;
nCpi = p.cpi*p.fs;

% fBand = .25e6;
% ds = 10;
winLen = 300;
nTaps = 20;

% p.fs = p.fs; % in order to remove
% cfs = (-11:2:11)/24; % frequency of channels
% profile on
for i = 1:floor(length(rxa)/nCpi)
  idx = (1:nCpi)+nCpi*(i-1);
  tCh = gpuArray(txa(idx));
  rCh = gpuArray(rxa(idx));
%   

% figure(56)
% hold off
% plot(log(smpl(real(abs(fftshift(fft(tCh)))),1000)))
% hold on
% plot(log(smpl(real(abs(fftshift(fft(rCh)))),1000)))


%   tChan = extractChannels(tCh,0,p.fs*ds,fBand);
%   rChan = extractChannels(rCh,0,p.fs*ds,fBand);
%   tCh = tChan(1:ds:end,:);
%   rCh = rChan(1:ds:end,:);
% alignFun(tCh,rCh);
%   [tCh,rCh] = alignFun(tCh,rCh);


% tCh = in;
% rCh = in2;

%   for j = 1
  [rxClean, tCh] = remDirect(winLen,nTaps,tCh,rCh);
  [r,v, rdi] = sledgeBatch(tCh,rxClean,p);
  
  plotRoutine(rdi,r,v,rCh(1:length(tCh)),rxClean,tCh)
%   
drawnow
%   end
  
  
end

save(num2str(p.fc/1e5),'rxa','txa')

% profile viewer 

% plot(abs(fftshift(fft(rx))))


























% alignFun(r(end-1e5:end),t(end-1e5:end));


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