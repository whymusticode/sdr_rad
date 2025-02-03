function [lag,TX,RX,MAX] = getLag(tx,rx,PLOT)


RX = fft(flipud(conj(rx)));
TX = fft(tx);
val = ifft(RX.*TX);
aval = fftshift(abs(val));

n = length(rx);
cpiLag = (1:n) - floor(n/2);

[MAX,loc] = max(aval);

lag = -cpiLag(loc);

if PLOT
%   figure(506)
plot(cpiLag,aval)
drawnow
end

