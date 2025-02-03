function out = extractChannels(in,cfs,bpFreq)

% in need to be column vector


% cfs = (-11:2:11)/24; % frequency of channels
% freqIn = 2.4e6; % sample freq in
% freqOut = .25e6; % band pass filter output
% bpFreq = .1 % ratio of sample frequency to band passed frequency



nCh = length(cfs); % number of channels (12)

N = length(in);

dConv = round(-cfs*N); % samples to shift to baseband


FILTER = fft(fir1(100,bpFreq/2)',N);

IN = fft(in);
out = zeros(N,nCh,'single','gpuArray');
for i = 1:nCh
  out(:,i) = ifft(circshift(IN,[dConv(i),0]).*FILTER);
end


% figure(23)
% plot(fftshift(abs(FILTER)))


% if dConv(i)>0
%     down = ifft([IN(end-dConv(i)+1:end);IN(1:end-dConv(i))].*FILTER);
%   else
%     down = ifft([IN(-dConv(i)+1:end);IN(1:-dConv(i))].*FILTER);
%   end
%   out(:,i) = down;





