

% clear

% 18500 watts (92.9 FM)
fc = 92.9e6; % Center frequency (Hz)
sampleRate  = 1e6;     % Samples per second
n = 2^10*256; % samples/frame (buffer size kindof)
N = 1e2; % number of frames
% n*N/sampleRate
% 3e8/fc/2   % bunny ears should be this size

sdrRx = comm.SDRRTLReceiver(...
  'CenterFrequency', fc, ...
  'EnableTunerAGC',  true, ...
  'SampleRate',      sampleRate, ...
  'SamplesPerFrame', n, ...%   'TunerGain', 20, ...
  'OutputDataType',  'double');
sdrRx2 = sdrRx;
sdrRx2.RadioAddress = '2';
% tic

DATA = gpuArray(zeros(n,N,2));


ds = 5; % downsample rate
KERN = gpuArray(fir1(1e2,1/ds));
lost = zeros(N,2);
% tic
for i = 1:N
  tic
%   [DATA(:,i,1), ~,lost(i,1)] = step(sdrRx);
%   toc
%   DATA(:,i) = DATA(:,i) - mean(DATA(:,i));
  [DATA(:,i,2), ~,lost(i,2)] = step(sdrRx2);
  toc
%   DATA(:,i,2) = DATA(:,i,2) - mean(DATA(:,i,2));
%   figure(3)
%   subplot(2,1,1)
%   plot(faf(DATA(:,i,1)))
%   font
%   subplot(2,1,2)
%   plot(faf(DATA(:,i,2)))
%   font
%   drawnow
%     step(hSpectrum, gather(DATA2(:,i)));
%     step(hSpectrum2, gather(DATA2(:,i)));
  %   DATA(:,i) = data;
  
%   step(hSpectrum,  step(hSDRrRx1));
end
% toc
if any(any(lost))
  disp('samples LOST !!!')
  disp(round([sum(lost1)/n/N,sum(lost2)/n/N]*100))
end
% profile viewer

% [val,lag] = xcorr(DATA(:,1),DATA2(:,1));
% plot(lag,abs(val))
% 
% return


return
dataOV = reshape(DATA,[numel(DATA)/2,2]);% OV for oversampled

data(:,1) = conv(dataOV(:,1),KERN,'same');
data(:,2) = conv(dataOV(:,2),KERN,'same');
dat = data(1:ds:end,:);


%% RD processing 
figure(1)
n2 = 1e4;
p.fs = sampleRate/ds; % sample frequency
p.fc = fc; % center frequency
p.rMax = 1e5; % maximum range
p.vMax = 600; % maximum doppler
tic
[r,v,rdi] = sledge(dat(1:n2,1),dat(1:n2,1),p);
toc
return

 figure(3)
  subplot(2,1,1)
  plot(faf(dataOV(:,1)))
  font
  subplot(2,1,2)
  plot(faf(dat(1:n2,2)))
  font

[val,lag] = xcorr(dat,dat2);
plot(lag,abs(val))
plot(1:length(dat),real(dat),1:length(dat2),real(dat2))
return

nnn = n;
tt = max(reshape(abs(data),[nnn,numel(data)/nnn]));
tt2 = max(reshape(abs(data2),[nnn,numel(data2)/nnn]));
% plot(1:length(tt),tt,1:length(tt),tt2)
corrcoef(tt,tt2)


C = 1e5;
d1 = data2(1:ds:C);
err = gpuArray(zeros(100,1));
I = 1:20:1e4;
for j = 1:length(I)
  i = I(j);
  d2 = data(i:ds:C+i-1);
  
  g = d1-d2;
  err(j) = real(g'*g);
end

figure(24)
% plot(1:C/ds,real(d1),1:C/ds,real(d2))


subplot(1,1,1)

% data = data(1:ds:end);
% data2 = data2(1:ds:end);
% nn = length(data);
plot(1:length(err),err)
font


% datLow = conv(DATA(:),KERN,'same');
% x = -100:100;
% plot(x,real(datLow(x+n)),x,imag(datLow(x+n)))
% plot(real(datLow(1:2*n)))


% plot(faf(datLow(1:2*n)))

% end
% toc
% 476 giga computations / sec    on GPU
% .6 giga computations / sec      on CPU

% 4.372 terraflops GPU theoretical
% 29.69 GFLOPS CPU theoretical






% kern = fir1(1e3,.1);
% tic
% datLow3 = conv(DATA(:),kern,'same');
%
%
% toc


% toc

% figure(9)
% i = randi(N-1);
% subplot(1,1,1)
% plot(-100:100,real([DATA(end-100:end,i);DATA(1:100,i+1)]))
% subplot(2,1,2)
% plot(fftshift(abs(fft(DATA(:,i)))))

% Release all System objects
% release(hSDRrRx);
% release(hSpectrum);

%% other shit
% timing on samples   
% subtract direct from look channel (regression model?)
% Compress look with direct (sledgehammer processing?)


