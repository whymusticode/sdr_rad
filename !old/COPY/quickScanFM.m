

% scans over FM frequencies and provides a spectrum 
clear
n = 1e4;
fs = 2.4e6;
fc = 87.9e6;
g = 20;
sdrRx = comm.SDRRTLReceiver(...
  'CenterFrequency', fc, ...
  'EnableTunerAGC',  false, ...% false enables 'TunerGain'
  'SampleRate',      fs, ...
  'SamplesPerFrame', n, ...
  'TunerGain',g, ...
  'OutputDataType',  'single');
sdrRx.RadioAddress = '0';

i = 0;
vals = [];
freqs = [];
while sdrRx.CenterFrequency < 108e6
i = i+1;

disp(round( sdrRx.CenterFrequency/1e6))

[DATA,~,~] = step(sdrRx);
vals(i,:) = smpl(abs(fftshift(fft(DATA))),10);

freqs(i,:) = smpl(((-floor(n/2):ceil(n/2)-1)/n*fs + ...
  sdrRx.CenterFrequency)/1e6,10);


release(sdrRx)

sdrRx.CenterFrequency = sdrRx.CenterFrequency + fs*.8;

end


figure(10)
plot(freqs',20*log10(vals'))
grid on

save('FMscan','freqs','vals')






