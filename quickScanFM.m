

% scans over FM frequencies and provides a spectrum 
% clear

n = 1e5;
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
DATA = DATA - mean(DATA);
vals(i,:) = smpl(abs(fftshift(fft(DATA))),10);

freqs(i,:) = smpl(((-floor(n/2):ceil(n/2)-1)/n*fs + ...
  sdrRx.CenterFrequency)/1e6,10);


release(sdrRx)

sdrRx.CenterFrequency = sdrRx.CenterFrequency + fs;

end


figure(13)
plot(freqs',20*log10(vals'))
grid on

save('FMscan','freqs','vals')


fr = freqs';fr = fr(:);
va = vals';va = va(:);

f = 3e2./dataM2(:,3);

for i = 1:length(f)
  idx = abs(f(i)-fr) < .12;
  thing(i) = sum(va(idx));
end
thresh =  1e5;
f(thing > thresh)

sqrt(val(thing > 1e4))


% dataM2(thing > 1e3,[1,2,6]) % location of emitters


%% reference:
% dataMat = [x,y,lam,pow,direc,zUp];












