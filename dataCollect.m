
% this is for setting gain values and storing off a certain amount of time
% at optimal gain values 

clear
% reset(sdrRx)
% reset(sdrTx)
% 87.9 - 107.9  %FM
% 389 + 6*rfChannel % DTV

fs = .25e6;
N = 375e3;

FC = 87.9:.2:107.9;
time = 1; % seconds of data to collect
buffer = 1; % seconds max offset between channels

goal = -3; % (db)
lower = -8; % (db)
glim = 40; % (db)

steps = 2;%round((time+buffer)*fs/N);
n = round(time*fs);


% data = zeros(n,2,length(FC),'int8');
for i = 1:length(FC)% find(FC==91.9)% 1:length(FC) % find(FC==91.9)%
  disp(i)
  fc = round(FC(i)*1e6);
  g = 10; % amplification at the beginning (db)
  
  clear sdrRx sdrTx
  sdrRx = comm.SDRRTLReceiver('CenterFrequency', fc, 'SampleRate',fs, ...
    'SamplesPerFrame',N,'OutputDataType',  'single'...
    ,'EnableTunerAGC',  false,'TunerGain',g);
  sdrTx = comm.SDRRTLReceiver('CenterFrequency', fc, 'SampleRate',fs,...
    'SamplesPerFrame',N,'RadioAddress','01','OutputDataType', ...
    'single','EnableTunerAGC',  false,'TunerGain',g);
  
  rx = []; tx = []; ct = 0;
  while ct < steps
    try
      tmp = step(sdrRx);
      tmp2 = step(sdrTx);
    catch
      rx = []; tx = []; ct = 0;
      tmp = step(sdrRx);
      tmp2 = step(sdrTx);
      disp('try statement worked')
    end
    
    val = 20*log10(max(abs(real([tmp;tmp2]))));
    disp([g,val])
    
    
    if val < lower && g < glim
      g = g + goal - val;
      g = double(g);
      release(sdrRx)
      release(sdrTx)
      sdrRx.TunerGain = g;
      sdrTx.TunerGain = g;
      rx = []; tx = []; ct = 0;
    elseif val == 0 && ct < 4
      g = g - 5; % roll off
      g = double(g);
      release(sdrRx)
      release(sdrTx)
      sdrRx.TunerGain = g;
      sdrTx.TunerGain = g;
      rx = []; tx = []; ct = 0;
    else % no problem || val < lower && g >= glim (tx is too weak/far away)
      ct = ct + 1;
      rx = [rx;tmp2];
      tx = [tx;tmp];
    end
    
  end
  
  
  %   playRx(vec2mat(rx,125*96).', fs)
  %   playRx(vec2mat(tx,125*96).', fs)
  [rxa,txa] = alignFun(rx,tx);
%     drawnow
%   if length(rxa) >= n
%     data(:,:,i) = [int8(txa(end-n+1:end)*2^8),int8(rxa(end-n+1:end)*2^8)];
%   end
  GAIN(i) = g;
end

store = [int8(txa(end-n+1:end)*2^8),int8(rxa(end-n+1:end)*2^8)];

save('959data','store')


%     if ct > steps % finally, enough data to store
%       [rxa,txa] = alignFun(rx,tx);
%       data(:,:,i) = [rxa(1:n),txa(1:n)];
%       continue
%     end

% clearvars -except data
% [rxa,txa,~,snr] = alignFun(data(:,1,1),data(:,2,1));


% save('ALLdata3','data','-v7.3')



% p.rMax = 1e5; % maximum range
% p.vMax = 600; % maximum doppler
% p.cpi = 2; % seconds to integrate


% return


% tic
% % rx = gather(rx);
% % tx = gather(tx);
% if check2 > check1 % use the sdr using more dynamic range as recieve
%   [txa,rxa] = alignFun(rx,tx);
% else
%   [rxa,txa] = alignFun(rx,tx);
% end
% toc

% alignFun(rxa(1:2e4),txa(1:2e4));
% rdi = LSsledge(txa(1:2e4),rxa(1:2e4),40,50);


% figure(461)
% imagesc(20*log10(abs(rdi')))
% colorbar

% return
% nCpi = p.cpi*p.fs;
%
% ds = 1;
% winLen = 300;
% nTaps = 20;

% p.fs = p.fs/ds; % in order to remove
% cfs = (-11:2:11)/24; % frequency of channels
% profile on
% for i = 1:floor(length(rxa)/nCpi)
%   idx = (1:nCpi)+nCpi*(i-1);
%   tCh = gpuArray(txa(idx));
%   rCh = gpuArray(rxa(idx));
%
%   tChan = extractChannels(in,0);
%   rChan = extractChannels(in2,0);
%   tCh = tChan(1:ds:end,:);
%   rCh = rChan(1:ds:end,:);

%   [tCh,rCh] = alignFun(tCh,rCh);


% tCh = in;
% rCh = in2;

%   for j = 1
%   [rxClean, tCh] = remDirect(winLen,nTaps,tCh,rCh);
%   [r,v, rdi] = sledgeBatch(tCh,rxClean,p);

%   plotRoutine(rdi,r,v,rCh(1:length(tCh)),rxClean,tCh)
%
% drawnow
%   end


% end
% profile viewer

% plot(abs(fftshift(fft(rx))))






% else
%   sdrRx = comm.SDRRTLReceiver('CenterFrequency', p.fc, 'SampleRate',p.fs, ...
%   'SamplesPerFrame',N,'OutputDataType',  'single'...
%   );
% sdrTx = comm.SDRRTLReceiver('CenterFrequency', p.fc, 'SampleRate',p.fs,...
%   'SamplesPerFrame',N,'RadioAddress','01','OutputDataType', ...
%   'single');
% end



















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