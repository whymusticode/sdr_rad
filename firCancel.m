function rxClean = firCancel(data,winLen,nTaps)
% finds least square fir filter applied to data(2:end) to cancel data(:,1)

rx = data(:,1);
n = size(rx,1);
nWins = floor(n/winLen);

%% now solve for FIR filters
H = [];
for i = 2:size(data,2)
  input = data(:,i);
  Hpre = repmat(input,[nTaps,1]);
  cols = floor(numel(Hpre)/(n+1));
  nT = round((cols-1)/2); % actual number of taps
  HH = reshape([input(end-nT+1:end);Hpre(1:(n+1)*cols-nT)],[n+1,cols]);
  HH = HH(1:end-1,:);
  H = cat(3,H,reshape(HH,[winLen,nWins,cols]));
end

H = permute(H,[1,3,2]);
Hp = conj(permute(H,[2,1,3]));
HH = mm3d(Hp,H);
RX2 = reshape(rx,[winLen,1,nWins]);

tmp = pagefun(@mldivide,HH,Hp);
ws = mm3d(tmp,RX2); % weights

direct = mm3d(H,ws); % the direct path

rxClean = rx - direct(:);

% figure(42)
% imagesc(abs(permute(ws,[3,1,2])))


% other ways of doing this:

% hhi = zeros(cols,cols,nWins,'single','gpuArray');
% for i = 1:nWins
%   hhi(:,:,i) = HH(:,:,i)^-1;%,Hp(:,:,i));
% end
% TX = mm3d(H,mm3d(mm3d(tmp,Hp),RX2));%arrayfun(@rdivide,HH,Hp)

% hhi = arrayfun(@inv,HH);
% tmp = mm3d(hhi,Hp);

% hhi = pagefun(@inv,HH);
% TX = mm3d(H,mm3d(mm3d(hhi,Hp),RX2));



%% EQ step (doesn't work)
% TXbl = fft(reshape(tx,[winLen,nWins]));
% TXfq = mean(abs(TXbl),2);
% RXfq = mean(abs(fft(reshape(rx,[winLen,nWins]))),2);
% txbl = ifft(TXbl.*RXfq./TXfq);
% tx = txbl(:);


% hold off;plot(log(fftshift(abs([TXfq,RXfq]))))
% hold off;plot(fftSmpl(rx));hold on;plot(fftSmpl(tx))
