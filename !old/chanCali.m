function out = chanCali(in,match)




n = length(match);

denom = fft(in(1:n));
% strength = median(abs(denom));

filter = ifft(fft(match)./(denom + 1*randiq((1:n)')));

n2 = length(in);
OUT = fft(in).*fft(filter,n2);
out = ifft(OUT);
