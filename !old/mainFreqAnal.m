function DATA = mainFreqAnal(sdrRx,N)
nSDRs = length(sdrRx);
n = sdrRx{1}.SamplesPerFrame;
fc = sdrRx{1}.CenterFrequency;
fs = sdrRx{1}.SampleRate;
freqs = ((-floor(n/2):ceil(n/2)-1)/n*fs + fc)/1e6;
DATA = zeros(n,nSDRs,N,'single');

figure(45)
for i = 1:N
  for j = 1:nSDRs
    [DATA(:,j,i), ~,lost(i,j)] = step(sdrRx{j});
    disp(length(unique(real(DATA(:,j,i)))))
    DATA(:,j,i) = DATA(:,j,i) - mean(DATA(:,j,i));
    tp = 20*log10(faf(DATA(:,j,i)));
    subplot(nSDRs,1,j)
    plot(freqs(1:100:end),smpl(tp,100));
    font
    drawnow
  end
end
if max(lost(:)) ~=0
  warning('lost samples')
end
% figure(46)

% plot(real(DATA(1:1e3,j,i)))
