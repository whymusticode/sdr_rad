

% scans over FM frequencies and provides a spectrum 
clear
n = 1e5;
fs = 2.4e6;
fc = 85e6;
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
while sdrRx.CenterFrequency < 110e6
i = i+1;

disp(round( sdrRx.CenterFrequency/1e6))

[DATA(:,i),~,~] = step(sdrRx);
vals(i,:) = smpl(abs(fftshift(fft(DATA(:,i)))),10);

freqs(i,:) = smpl(((-floor(n/2):ceil(n/2)-1)/n*fs + ...
  sdrRx.CenterFrequency)/1e6,10);


release(sdrRx)

sdrRx.CenterFrequency = sdrRx.CenterFrequency + fs;

end


figure(11)
plot(freqs',20*log10(vals'))
grid on

save('FMscan','freqs','vals')






