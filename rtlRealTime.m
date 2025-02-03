% 
% clear
% 
% nTaps = 6; % number of taps in cancellation filter
% winLen = 100; % length of cancellation window
% p.rMax = 2e5; % maximum range
% p.vMax = 600; % maximum doppler
% p.cpi = 1; % seconds to integrate
% fc = 98e6; % Center frequency (Hz)
% fs  = 3.2e6; %(hz) 225e3 < fs <= 300e3 or 900e3 < fs <= 3200e3
% n = 375e3; % samples/frame  [maximum is 375e3]
% g = 15; % amplification at the beginning (db?)
% buffer = 4*fs; % maximum number of samples channels might be off
% Nsdr = 3; % number of SDRs connected
% 
% 
% tRun = 80; % code taking too long to execute at > 80, 70 runs well 
% % ('EnableBurstMode' has no effect on this)
% % variable in memory where output is stored must be less than 1 GB or
% % stepping SDRs starts taking too long
% 
% Texp = n/fs
% Nct = round(tRun/Texp);
% % data = zeros(n,Nct,Nsdr,'int16','gpuArray');
% lost = zeros(Nct,Nsdr);
% 
% %% make SDR objects
% % for k = 1:length(fc)
% for i = 1:Nsdr
%   sdrRx{i} = comm.SDRRTLReceiver(...
%     'CenterFrequency', fc, ...
%     'EnableTunerAGC',  true, ...% 'TunerGain',g, ...
%     'SampleRate',      fs, ...
%     'SamplesPerFrame', n, ...
%     'RadioAddress',num2str(i-1),...
%     'OutputDataType',  'int16');
% %   'EnableBurstMode', true,...
% end
% 
% for ct = 1:Nct
% % ct = 0;
% % while true
% %   ct = ct+1
%   tic
%   for j = 1:Nsdr
%     [data,~,lost(ct,j)] = step(sdrRx{j});
% %     disp(max(real(data(:,ct,j))))
%   end
%   toc
% end
% %   save(num2str(fc(k)),'-v7.3')
% % end
% % 20*log10(double(max(real(data(:))))/128)
% CH1 = data(:,:,2);
% CH2 = data(:,:,3);
% 
% % offset = 950934;
% 
% 
% 
% 
% testLen = 1e6;
% 
% i = 0
% 
% for i = 1:100
% ch1 = gS(CH1((1:testLen)+offset + i*testLen));
% ch2 = gS(CH2((1:testLen) + i*testLen));
% 
% [lag,TX,RX,MAX] = getLag(ch1(:),ch2(:),1);
% val(i) = lag;
% val2(i) = MAX;
% % offset = offset - lag;
% title(i)
% pause(.02)
% end
% % drawnow
% 
% 
% plot(val2,'.')
% 
% magdb = 20*log10(abs(fftshift(fft(CH1(1:1e6)))));
% ans = smpl(magdb,100);
% freqs = (freqBins(length(ans))*fs + fc(2))/1e6;
% plot(freqs,ans)



load('973_oneChannel.mat')
dat = reshape(data,size(data,1)*size(data,2),size(data,3));
for i = 1:2
  offsets(i) = getLag(gS(dat(1:1e7,4-i)),gS(dat(1:1e7,3-i)),0);
end
% offsets = [955001,950934];
datAligned = gS([dat(sum(offsets)+1:end,1),...
  dat(offsets(1)+1:end-offsets(2),2),dat(1:end-sum(offsets),3)]);

% getLag(datAligned(1:1e7,2),datAligned(1:1e7,3),1);
% xlim([-1,1]*100)



winLen = 100;
nTaps = 20;
dat2 = datAligned((1:1e5)+2e6,[3,2]);

rxClean = firCancel(dat2,winLen,nTaps);

20*log10(real((rxClean'*rxClean)/(dat2(:,1)'*dat2(:,1))))








