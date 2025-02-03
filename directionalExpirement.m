%#ok<*UNRCH>

% scan all FM radio channels ot check for how strong the singal is
clear

% release(sdr1)
% release(sdr2)

fs = .25e6;
N = 375e3;
FC = 91.9; 
g = 25; % amplification of bunny ear antenna
g2 = 25; % gain of log periodic antenna 

fc = round(FC*1e6);
sdr1 = comm.SDRRTLReceiver('CenterFrequency', fc, 'SampleRate',fs, ...
  'SamplesPerFrame',N,'OutputDataType',  'single'...
  ,'EnableTunerAGC',  false,'TunerGain',g);
sdr2 = comm.SDRRTLReceiver('CenterFrequency', fc, 'SampleRate',fs,...
  'SamplesPerFrame',N,'RadioAddress','01','OutputDataType', ...
  'single','EnableTunerAGC',  false,'TunerGain',g2);

data = [];
while 1
  
  tmp = step(sdr1);
  tmp2 = step(sdr2);
  data = cat(1,data,[tmp,tmp2]);
  
  strength = 20*log10(max(abs(real([tmp,tmp2]))));
  disp(strength)
end


[rx,tx,smallLag,snr] = alignFun(data(:,2),data(:,1)); 
data = [tx,rx];
data = int8(data*2^7);
save('919BunnyEarDirect','data')