function [rx,tx,smallLag,snr] = alignFun(rx,tx)
n = length(rx);

maxLen = 2e7;
len2use = min(maxLen,n);
cpiLag = (1:len2use) - floor(len2use/2);

ccf = abs(fftshift(ifft(fft(rx(1:len2use)).*fft(conj(flipud(tx(1:len2use)))))));

[val,loc] = max(ccf);
snr = val/median(ccf);
smallLag = cpiLag(loc);
if smallLag > 0; rx = rx(smallLag+1:end);  tx = tx(1:end-smallLag);
else; rx = rx(1:end+smallLag);  tx = tx(-smallLag+1:end); end

figure(409)
plot(smpl(cpiLag,1e3),smpl(ccf,1e3))
title('alignment')

% do math
% if isempty(varargin)

% else
% buf = 1e8; % buffer, not always necessary
% integ = maxN*2;
% ccf = abs(fftshift(ifft(fft(rx((1:integ)+buf)).*...
%   fft(conj(flipud(tx((1:integ)+buf)))))));
% end

% plot(cpiLag,ccf)
%