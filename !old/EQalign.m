

% clear
% reset(gpuDevice)
% load('highSampleRateTimeAlignedData.mat')


intTime = 2;
nCpi = intTime*p.fs;
cfs = ((2:2:22)-12)/24;
chBW = .3e6;



RX = fftshift(fft(rxa(1:nCpi)));
TX = fftshift(fft(txa(1:nCpi)));
% RX = fft(rxa(1:1e6));
% TX = fft(txa(1:1e6));

x = smpl(p.fc + p.fs*freqBins(nCpi),100)/1e6;
figure(20)
subplot(3,1,1)
plot(x,log(smpl(RX,100)),x,log(smpl(TX,100)))
grid on

out = extractChannels(rxa(1:nCpi),cfs,p.fs,chBW);
subplot(3,1,2)
hold off
for i = 1:length(cfs)
  small = out(1:8:end,i);
  y = log(abs(fftshift(fft(small))));
  x = freqBins(length(y))*chBW + cfs(i)*p.fs + p.fc;
  plot(smpl(x,100)/1e6,smpl(y,100))
  hold on
end
grid on


out2 = extractChannels(txa(1:nCpi),cfs,p.fs,chBW);
subplot(3,1,3)
hold off
for i = 1:length(cfs)
  small = out2(1:8:end,i);
  y = log(abs(fftshift(fft(small))));
  x = freqBins(length(y))*chBW + cfs(i)*p.fs+ p.fc;
  plot(smpl(x,100)/1e6,smpl(y,100))
  hold on
end
grid on

figure(4)
tic
for i = 1:11
rx = out(1:8:end,i);
tx = out2(1:8:end,i);
rxClean = removeDirectVect(200,30,tx,rx);
disp(10*log10(real(rxClean'*rxClean/(rx'*rx))))

plot([fftSmpl(rx);fftSmpl(tx);fftSmpl(rxClean)]')
font
legend('rx','tx','rxClean')
pause(1)

drawnow
end
toc


