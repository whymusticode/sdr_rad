  
%   subplot(2,1,1)
%   DATA(:,j,1) = DATA(:,j,1) - mean(DATA(:,j,1));
%   tp = 20*log10(faf(DATA(:,j,1)));
%   plot(smpl(gather(tp)));
%   font
%   drawnow
%   subplot(2,1,2)
%   DATA(:,j,2) = DATA(:,j,2) - mean(DATA(:,j,2));
%   tp = 20*log10(faf(DATA(:,j,2)));
%   plot(smpl(gather(tp)));
%   font
%   drawnow
%   
%   disp([length(unique(DATA(:,j,1))),length(unique(DATA(:,j,2)))])

% plot(1:n,real(DATA(:,1)),1:n,real(DATA(:,2)))

% disp([length(unique(real(DATA(:,1)))),length(unique(real(DATA(:,2))))])



return
n = floor(length(DATA)/650)*650;
DATA = DATA(1:n);

data = conv(DATA,gather(KERN),'same');
[z,zbb] = fmDemod(data);
y = fmMod(z,zbb);

plot(1:n,real(DATA(:,1)),1:n,real(DATA(:,2)))
xlim([1e5,1e5+1e4])
% plot(log(faf(DATA)))
[val,lag] = xcorr(data,y);
plot(lag,abs(val))


return

%% demodulation remodulation idea

% fDev = 75000;
% filterTimeConst = 7.5e-5;
% fmDemod = comm.FMBroadcastDemodulator(...
%   'SampleRate', fs, 'FrequencyDeviation', ...
%   fDev, 'FilterTimeConstant', filterTimeConst,...
%   'Stereo', true,'RBDS',true,'AudioSampleRate', audioSampleRate);
% fmMod = comm.FMBroadcastModulator('SampleRate',fs,...
%   'FrequencyDeviation',fDev,'Stereo', true,'RBDS',true, ...
%   'FilterTimeConstant', filterTimeConst,'AudioSampleRate',audioSampleRate);

% z = fmdemod(DATA,fc,fs,fDev); % Demodulate both channels.
% y = fmmod(x,fc,fs,fDev);




% groupLen = 104;
% sps = 10;
% groupsPerFrame = 19;
% rbdsFrameLen = groupLen*sps*groupsPerFrame;
% afrRate = 40*1187.5;
% rbdsRate = 1187.5*sps;
%
% afr = dsp.AudioFileReader('rbds_capture_47500.wav','SamplesPerFrame',...
%   rbdsFrameLen*afrRate/rbdsRate);
% rbds = comm.RBDSWaveformGenerator('GroupsPerFrame',groupsPerFrame,...
%   'SamplesPerSymbol',sps);
%
% fmMod = comm.FMBroadcastModulator('AudioSampleRate',afr.SampleRate,...
%   'SampleRate',fs,...
%     'Stereo',true,'RBDS',true,'RBDSSamplesPerSymbol',sps);
% fmDemod = comm.FMBroadcastDemodulator('SampleRate',outRate,...
%     'Stereo',true,'RBDS',true,'PlaySound',true);
% scope = dsp.TimeScope('SampleRate',outRate,'YLimits',10^-2*[-1 1]);
%
% for idx = 1:7
%     input = afr();
%     rbdsWave = rbds();
%     yFM = fmMod([input input], rbdsWave);
%     rcv = awgn(yFM, 40);
%     [audioRcv, rbdsRcv] = fmDemod(rcv);
% yFM2 = fmMod(audioRcv(1:length(input),:), rbdsRcv);
%     scope(rbdsRcv);
% end
% hold off
% plot(input)
% hold on
% plot(audioRcv)

% fs = 100;              % Sample rate (Hz)
% ts = 1/fs;             % Sample period (s)
% fd = 25;               % Frequency deviation (Hz)
% t = (0:ts:0.5-ts)';
% x = sin(2*pi*4*t);
% y = fmMod(x);
% z = fmDemod(y);
% plot(t,x,t,z)


% clear
% fmrxCust(10, 92.9)


% load('coherentMinus3.mat')


rcv = gather(DATA(1:floor(end/(650*96))*650*96,1,1));
[audioRcv, rbdsRcv] = fmDemod(rcv');
rcvEst = fmMod(rbdsRcv);



return
DATA = gpuArray(zeros(n,1));

KERN = gpuArray(fir1(1e2,1/ds));
lost = zeros(N,1);
for i = 1:N
  [DATA, ~,lost(i)] = step(sdrRx);
  DATA = DATA - mean(DATA);
end
% len = floor(length(DATA)/650)*650;
% rcvDmd = fmDemod(DATA(1:len));
% rcvEst = fmMod(rcvDmd);
data = conv(DATA,gather(KERN),'same');
dat = data(1:ds:end,:);

% can we reconstruct the signal?

fcb = fDev*2;
modu = exp(1i*2*pi*(1:n)*fcb/fs)';

x = real(data.*modu);
z = fmdemod(x,fcb,fs,fDev);
y = fmmod(z,fcb,fs,fDev);


toc
hold off
plot(real(x./modu))
hold on
plot(real(y./modu))






return
plot(real(DATA(:,1)))
plot(10*log10(faf(DATA(:,1))))


[audioRcv, rbdsRcv] = fmDemod(gather(dat(1:floor(length(dat)/130)*130)))

[audioRcv, rbdsRcv] = fmDemod(rcv);
rcvEst = fmMod(rbdsRcv);

plot(real(dat))
plot(10*log10(faf(dat)))
% toc
if any(any(lost))
  disp('samples LOST !!!')
  disp(round([sum(lost1)/n/N,sum(lost2)/n/N]*100))
end
return

% sdrRx2 = sdrRx;
% sdrRx2.RadioAddress = '2';
% tic


%   toc
%   DATA(:,i) = DATA(:,i) - mean(DATA(:,i));
%   [DATA(:,i,2), ~,lost(i,2)] = step(sdrRx2);
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

% profile viewer


%% subtract the channels
[val,lag] = xcorr(DATA(:,3,1),DATA(:,7,2));
[~,init] = max(abs(val));
plot(lag,abs(val))
locInit = lag(init);

% [val,lag] = xcorr(d1,d2);
% [~,locInit] = max(abs(val));
% plot(lag,abs(val))


d2 = reshape(DATA(:,6:10,2),[5*size(DATA,1),1]);
d1 = reshape(DATA(:,2:6,1),[5*size(DATA,1),1]);
d1sm = conv(d1,KERN,'same');
d2sm = conv(d2,KERN,'same');

idx = round(length(d1sm)/2) + (0:1e4);




[val,lag] = xcorr(d1sm(idx),d2sm(idx+locInit));
figure(1)
hold off
plot(lag,abs(val))
figure(1)
hold off
plot(real(d1sm(idx)))
hold on
plot(real(d2sm(idx-locInit)))

% d1sm(round(end/2):

mean(abs(d1sm(idx)).^2)

ID = 1:length(d2);
TRY = -40:20;
for i = 1:length(TRY)
  
  sig2 = interp1(ID,d2sm,idx-locInit-TRY(i));
  VAL(i) = mean(abs(d1sm(idx) - conj(sig2')).^2);
end
figure(2)
plot(TRY,VAL)
font
% plot(lag,abs(val))

plot(real(df))
hold on
plot(real(df2))

df = d1sm(locInit:ds:end);
df2 = d2sm(1:ds:end-locInit);
clear DATA d1 d2 d1sm d2sm

p.fs = fs/ds; % sample frequency
p.fc = fc; % center frequency
p.rMax = 1e5; % maximum range
p.vMax = 300; % maximum doppler
[r,v,rdi] = sledge(df(1:end/2),df2(1:end/2),p);





%
% return

DATA(:,:,1) = DATA1(:,:,1);
return
% dataOV = reshape(DATA,[numel(DATA)/2,2]);% OV for oversampled
clear DATA
data(:,1) = conv(dataOV(:,1),KERN,'same');
data(:,2) = conv(dataOV(:,2),KERN,'same');
dat = data(1:ds:end,:);


%% RD processing
figure(1)
n2 = 1e4;
p.fs = fs/ds; % sample frequency
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


%% in recording loop to determine bitdepth
%   length(unique(real(DATA)))
%
%   plot(10*log10(faf(DATA)))
%   ylim([10,50])
%   drawnow
%   toc