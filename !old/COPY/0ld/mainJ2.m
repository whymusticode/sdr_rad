

% reset(gpuDevice(1)) % this will free up memory
clear
% fmrxCust(10, 92.9)
% 18500 watts (92.9 FM)
fc = 92.9e6; % Center frequency (Hz)
fs  = 2.6e6; % (hz) % success at 2.6 and tunergain at 100 (max 3.2e6)
n = 375e3; % samples/frame  [maximum is 375e3]
N = 10; % number of CPIs
% n*N/sampleRate
% 3e8/fc/4*39.4% bunny ears should be a little less than this size (inches)

sdrRx = comm.SDRRTLReceiver(...
  'CenterFrequency', fc, ...
  'EnableTunerAGC',  true, ...
  'SampleRate',      fs, ...
  'SamplesPerFrame', n, ...%   'TunerGain',100, ...
  'RadioAddress','01',...
  'OutputDataType',  'single');

% mainFreqAnal(sdrRx,N)
% return

ds = round(fs/256e3/1.5); % downsample rate
lost = zeros(2,N); % lost data
DATA = zeros(n,N,'single','gpuArray');


% figure(45)
% subplot(1,1,1)

for i = 1:N
    [DATA(:,i), ~,lost(i)] = step(sdrRx);
%     disp(length(unique(real(DATA(:,i,j)))))
%     DATA(:,i,j) = DATA(:,i,j) - mean(DATA(:,i,j));
%     tp = 20*log10(faf(DATA(:,i,j)));
%     subplot(2,1,j)
%     plot(smpl(gather(tp)));
%     font
%     drawnow
end
if max(lost(:)) ~=0
  error('lost samples')
end
return
figure(46)
% plot(real(DATA(1:1e3,j,i)))
data = DATA(:,1,1);
data2 = DATA(:,1,2);
[val,lag] = xcorr(data,data2);
plot(lag,abs(val))



return



for i = 1:2 % sample on both SDRs
  for j = 1:N
    [DATA(:,i,j), ~,lost(i,j)] = step(sdrRx{i});
  end
end

plot(real(DATA(:,1,1)))










