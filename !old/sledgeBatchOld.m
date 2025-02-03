function [range,velocity,image] = sledgeBatch(tx,rx,p)
% performs 2-D Cross Correlation Function in range and doppler for
% arbitrary tansmit waveform. Is very fast
%
% inputs: 
%   tx - IQ signal transmitted
%   rx - IQ signal recieved length(rx) == length(tx)
%   p - parameters including:
%     p.fs      sample frequency (Hz)
%     p.fc      center frequency (Hz)
%     p.rMax    maximum range to create image (m)
%
% outputs:
%   range - (m)
%   velocity - (m/s)
%   image - complex range doppler image, columns are range rows are 
%           velocity, like god intended 
% 
% matches up with imagesc format:  
% [x,y,z] = sledgeBatch(tx,rx,p); imagesc(x,y,z)


c = 3e8; % m/s speed of light
rpp = c/p.fs; % range per pixel
range = rpp:rpp:p.rMax; % range bins
nR = length(range); % batch length
nB = floor(length(rx)/nR); % number of batches
n = nR*nB; % truncated integration 
dpp = c/n*p.fs/p.fc; % doppler per pixel
velocity = (-floor(nB/2):ceil(nB/2)-1)*dpp; % velocity bins 

% circshift hack makes it cover the full range, this noncoherently sums 2
% different noise floors

SR = fft(reshape(rx(1:n),[nR,nB]));
SREF = fft(flipud(conj(reshape(tx(1:n),[nR,nB]))));
corFun1d = ifft(SR.*SREF); % range  + circshift(SREF,[0,1]))
image = fftshift(fft(corFun1d,[],2),2)'; % doppler 

% corFun1d = ifft(fft(rx(1:n)).*fft(flipud(conj(tx(1:n))))); % range
% batches = reshape(corFun1d,[nR,nB]);
% image = fftshift(fft(batches,[],2),2)'; % doppler 


% corFun1d2 = ifft(fft(sr).*fft(circshift(flipud(sref),[0,1]))); % range

% a = randn(4);
% fft(fft(a),[],2)
% fft2(a)








