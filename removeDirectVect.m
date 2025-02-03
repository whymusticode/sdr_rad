function rxClean = removeDirectVect(winLen,nTaps,tx,rx) %RLSest,LMSest,

% n = (length(rx)/winLen)



n = numel(tx);
nWins = floor(n/winLen);


%% EQ step (doesn't work)
% TXbl = fft(reshape(tx,[winLen,nWins]));
% TXfq = mean(abs(TXbl),2);
% RXfq = mean(abs(fft(reshape(rx,[winLen,nWins]))),2);
% txbl = ifft(TXbl.*RXfq./TXfq);
% tx = txbl(:);


% hold off;plot(log(fftshift(abs([TXfq,RXfq]))))
% hold off;plot(fftSmpl(rx));hold on;plot(fftSmpl(tx))


%% now solve for FIR filters
input = tx;
Hpre = repmat(input,[nTaps,1]);
cols = floor(numel(Hpre)/(n+1));
nT = round((cols-1)/2); % actual number of taps 
HH = reshape([input(end-nT+1:end);Hpre(1:(n+1)*cols-nT)],[n+1,cols]);
HH = HH(1:end-1,:);
H = reshape(HH,[winLen,nWins,cols]);
H = permute(H,[1,3,2]);
Hp = conj(permute(H,[2,1,3]));
HH = mm3d(Hp,H);
RX2 = reshape(rx,[winLen,1,nWins]);

tmp = pagefun(@mldivide,HH,Hp);
ws = mm3d(tmp,RX2); % weights
% ws = repmat(ws(:,1,1),[1,1,size(ws,3)]);
% figure(42)
% imagesc(abs(permute(ws,[3,1,2])))

direct = mm3d(H,ws); % the direct path 

rxClean = rx - direct(:);




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


