function [range,velocity,image] = sledgeBatch(tx,rx,p)
% performs 2-D Cross Correlation Function in range and doppler for
% arbitrary tansmit waveform. Caution, this has doppler ambiguities
%
% inputs:
%   tx - IQ signal transmitted
%   rx - IQ signal recieved length(rx) == length(tx)
%   p - parameters including:
%     p.fs      sample frequency (Hz)
%     p.fc      center frequency (Hz)
%     p.rMax    maximum range to create image (m)
%     p.vMax    max velocity +- to create image (m/s)
%
% outputs:
%   range - (m)
%   velocity - (m/s)
%   image - complex range doppler image, [velocity x range]
%
% matches up with imagesc format:
% [r,v,rdi] = sledgeBatch(tx,rx,p); imagesc(r,v,rdi)


% create bins
c = 3e8; % m/s speed of light
n0 = length(tx); % length of inputs
dpp = c/n0*p.fs/p.fc; % doppler per pixel approximation
nD = round(2*p.vMax/dpp); % number of doppler bins
nB = floor(n0/nD); % batch length
n = nD*nB; % new integration length
dpp = c/n*p.fs/p.fc; % doppler per pixel exact
rpp = c/p.fs; % range per pixel
range = rpp:rpp:p.rMax; % range bins
nR = length(range); % number of range bins
velocity = (-floor(nD/2):ceil(nD/2)-1)*-dpp; % velocity bins
nBR = nB + nR; % length of convolution




% index 
idx = (1:nBR)' + (0:(nD-1))*nB;
idxidx = idx > n0;
idx(idxidx) = 1;
rxBatch = rx(idx);
rxBatch(idxidx) = 0;

% processing
TX = fft( flipud(conj( reshape(tx(1:n),[nB,nD]) )) ,nBR);
RX = fft(rxBatch); 
if 0
  corFun1d = ifft(RX.*TX.*fftshift(hamming(nBR)));
  image = fftshift(fft(corFun1d(nB+1:end,:).*hamming(nD)',[],2),2)';
else
  corFun1d = ifft(RX.*TX);
  image = fftshift(fft(corFun1d(nB+1:end,:),[],2),2)';
end


