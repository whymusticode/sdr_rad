
reset(gpuDevice)
clear

p.fs = .3e6;
N = 375e3; % samples of data pulled on every step

p.rMax = 1e5; % maximum range
p.vMax = 600; % maximum doppler
p.cpi = 2; % seconds to integrate

g = 40;
gLim = 50;


winLen = 200;
nTaps = 30;

time = 30; % seconds of data to collect

% FC = 87.9:.2:107.9; % all broadcast center frequencies
FC = [89.7,92.5,92.9,93.7,94.5,96.9,98.3,98.5,99.5,100.7,...
  102.5,103.3,104.1,105.7,106.7,107.3];

steps = round(time*p.fs/N);
nCpi = p.cpi*p.fs;

cancelPerf = zeros(size(FC));
tic
for j = 1:length(FC)
  p.fc = round(FC(j)*1e6);
  filename = [num2str(p.fc/1e5),'DataCollect7_29'];
  
  disp(num2str(FC(j)))
  
  sdrRx = comm.SDRRTLReceiver('CenterFrequency', p.fc, 'SampleRate',...
    p.fs, 'SamplesPerFrame',N,'OutputDataType',  'single',...
    'EnableTunerAGC',  false,'TunerGain',g);
  sdrTx = comm.SDRRTLReceiver('CenterFrequency', p.fc, 'SampleRate',...
    p.fs,'SamplesPerFrame',N,'RadioAddress','01','OutputDataType', ...
    'single','EnableTunerAGC',  false,'TunerGain',g);
%   sdrRx = comm.SDRRTLReceiver('CenterFrequency', p.fc, 'SampleRate',...
%     p.fs, 'SamplesPerFrame',N,'OutputDataType',  'single',...
%     'EnableTunerAGC',  true);
%   sdrTx = comm.SDRRTLReceiver('CenterFrequency', p.fc, 'SampleRate',...
%     p.fs,'SamplesPerFrame',N,'RadioAddress','01','OutputDataType', ...
%     'single','EnableTunerAGC', true);
  
  
  rx = zeros(steps*N,1);
  tx = rx;
  lost = zeros(steps,2);
  for i = 1:steps
    [rx((1:N)+N*(i-1)), sdrRx,lost(i,1),check1] = myStep(sdrRx);
    
%     if check1 < 100 % my version of automatic gain control
%       sdrRx.TunerGain = sdrRx.TunerGain + 2;
%     elseif check1 > 220
%       sdrRx.TunerGain = sdrRx.TunerGain - 2;
%     end
%     if sdrRx.TunerGain > gLim;sdrRx.TunerGain = gLim; end
%     
    [tx((1:N)+N*(i-1)), sdrTx,lost(i,2),check2] = myStep(sdrTx);
%     
%     if check2 < 100 % my version of automatic gain control
%       sdrTx.TunerGain = sdrTx.TunerGain + 2;
%     elseif check2 > 220
%       sdrTx.TunerGain = sdrTx.TunerGain - 2;
%     end
%     if sdrTx.TunerGain > gLim;sdrTx.TunerGain = gLim; end
  end
  %   disp(['lost this many samples: ',num2str(sum(lost(:)))])
  data.clock = clock;
  
  %   if check2 > check1 % use the sdr using more dynamic range as recieve
  [rxa,txa] = alignFun(rx,tx);
  %   else
  %     [rxa,txa] = alignFun(rx,tx);
  %   end
  
  rdi = [];
  for i = 1:floor(length(rxa)/nCpi)
    idx = (1:nCpi)+nCpi*(i-1);
    tCh = gpuArray(txa(idx));
    rCh = gpuArray(rxa(idx));
    
    
    %     [tCh,rCh] = alignFun(tCh,rCh);
    [rxClean, tCh] = remDirect(winLen,nTaps,tCh,rCh);
    [r,v, tmp] = sledgeBatch(tCh,rxClean,p);
    if isempty(rdi)
      rdi = tmp;
    else
      rdi = cat(3,rdi,tmp);
    end
    plotRoutine(rdi(:,:,i),r,v,rCh(1:length(tCh)),rxClean,tCh);
    drawnow
    
  end
  clear data
  data.clock = clock;
  data.rx = rxa;
  data.tx = txa;
  data.rdi = rdi;
  data.r = r;
  data.v = v;
  data.params = p;
  data.cancelPerf = 10*log10(abs(rxClean'*rxClean/(rCh'*rCh)));
  
    cancelPerf(j) = gather(10*log10(abs(rxClean'*rxClean/(rCh'*rCh))));
%   save(filename,'data')
  
  release(sdrRx)
  release(sdrTx)
end

toc

















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