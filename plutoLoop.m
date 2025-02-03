

clear

myRadios =  findPlutoRadio;


fc = 92.5e6; % Center frequency (Hz) [88.5]
fs  = .3e6; %(hz) 225e3 < fs <= 300e3 or 900e3 < fs <= 3200e3
n = 9600; % samples/frame  [maximum is 375e3] % 9600 for 300khz playback
g = 35; % amplification at the beginning (db)
secs = 3; % seconds of data to record 

steps = round(secs*fs/n);
% steps = 1;

sdr{1} = sdrrx('Pluto', ...
  'CenterFrequency', fc, ...
  'GainSource', 'Manual', ...
  'Gain', g, ...
  'BasebandSampleRate', fs,...
  'SamplesPerFrame',n,...
  'OutputDataType','single',...
  'RadioID', myRadios(1).RadioID);
sdr{2} = sdrrx('Pluto', ...
  'CenterFrequency', fc, ...
  'GainSource', 'Manual', ...
  'Gain', g, ...
  'BasebandSampleRate', fs,...
  'SamplesPerFrame',n,...
  'OutputDataType','single',...
  'RadioID', myRadios(2).RadioID);
release(sdr{1})
release(sdr{2})

data = zeros(n,steps,2);
for i = 1:steps
  [data(:,i,1),sdr{1}] = myStepFun(sdr{1});
  [data(:,i,2),sdr{2}] = myStepFun(sdr{2});
end


figure(67)
subplot(2,1,1)
plotFreq(data(:,1,1),fc,fs)
subplot(2,1,2)
plotFreq(data(:,1,2),fc,fs)


% playRx(data(:,:,1), fs)
% playRx(data(:,:,2), fs)

tx = data(:,:,1);
rx = data(:,:,2);
[lag,TX,RX,MAX] = getLag(tx(:),rx(:),1);

% rxClean = removeDirectVect(200,40,tx(:),rx(:));
return


rx = data(:,:,1);
tx = data(:,:,2);
rx = rx(:);
tx = tx(:);

[lag,TX,RX,MAX] = getLag(tx,rx,1);

p.rMax = 100e3;
p.vMax = 10e4;
p.fs = fs;
p.fc = fc;
[range,velocity,image] = ...
  sledgeBatch(tx(end-fs+1:end),rx(end-fs+1:end),p);

figure(402)
imagesc(range/1e3,velocity,log(abs(image)))

figure(405)  
subplot(2,1,1)
plot(fftshift(smpl(log(abs(fft(rx))),1e3)))
subplot(2,1,2)
plot(fftshift(smpl(log(abs(fft(tx))),1e3)))











%% end used code


% for j = 1:length(sdr)
%   release(sdr{j})
% end

% save params
% delete 'data\*.mat'


% addpath('data')





