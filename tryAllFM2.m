
% scan all FM radio channels ot check for how strong the singal is
clear


fs = .25e6;
N = 375e3;
FC = 87.9:.2:107.9;
g = 20; % amplification (db)

for i = 1:length(FC)% find(FC==91.9)% 1:length(FC) % find(FC==91.9)%
  disp(i)
  fc = round(FC(i)*1e6);
  
%   release(sdrRx)
%   release(sdrTx)
  sdrRx = comm.SDRRTLReceiver('CenterFrequency', fc, 'SampleRate',fs, ...
    'SamplesPerFrame',N,'OutputDataType',  'single'...
    ,'EnableTunerAGC',  false,'TunerGain',g);
  sdrTx = comm.SDRRTLReceiver('CenterFrequency', fc, 'SampleRate',fs,...
    'SamplesPerFrame',N,'RadioAddress','01','OutputDataType', ...
    'single','EnableTunerAGC',  false,'TunerGain',g);
  
  
  tmp = step(sdrRx);
  tmp2 = step(sdrTx);
  
  
  strength(i,:) = 20*log10(max(abs(real([tmp,tmp2]))));
end

save('channelStrength','strength','FC')