clear
fc = 88.5e6; % Center frequency (Hz)
fs  = .3e6; % (hz) % success at 2.6 and tunergain at 100 (max 3.2e6)
% can not be 3e5< fs <=9e5
n = 375e3; % samples/frame  [maximum is 375e3]
N = 300; % number of CPIs
% n*N/sampleRate
% 3e8/fc/4*39.4% bunny ears should be a little less than this size (inches)
g = 20;

sdrRx = comm.SDRRTLReceiver(...
  'CenterFrequency', fc, ...
  'EnableTunerAGC',  false, ...% set to false to manually tune
  'SampleRate',      fs, ...
  'SamplesPerFrame', n, ...%   'TunerGain',100, ...
  'TunerGain',g, ...
  'RadioAddress','01',...
  'OutputDataType',  'single');


% sdrRx.TunerGain;
i = 0;
while true
  i = i + 1;
  [DATA, ~,lost] = step(sdrRx);
  if lost ~=0;warning('lost samples in auxiliary');end
  check = length(unique(real(DATA)));
  disp(check)
  
  if check < 200 % my version of automatic gain control
    sdrRx.TunerGain = sdrRx.TunerGain + 3;
  elseif check > 250
    sdrRx.TunerGain = sdrRx.TunerGain - 3;
  end
  
  save(['aux',num2str(mod(i,3))],'DATA')
  
  tp = 20*log10(faf(DATA));
  plot(smpl(tp,100))
  font
  drawnow
  
end





