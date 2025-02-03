function [rxClean, tx] = remDirect(winLen,nTaps,tx,rx) 

n = floor(length(rx)/winLen)*winLen;
rx = rx(1:n);
tx = tx(1:n);



input = tx;
Hpre = repmat(input,[nTaps,1]);
cols = floor(numel(Hpre)/(n+1));
nT = round((cols-1)/2); % actual number of taps 
HH = reshape([input(end-nT+1:end);Hpre(1:(n+1)*cols-nT)],[n+1,cols]);
HH = HH(1:end-1,:);
nWins = floor(n/winLen);
H = reshape(HH,[winLen,nWins,cols]);
H = permute(H,[1,3,2]);
Hp = conj(permute(H,[2,1,3]));
HH = mm3d(Hp,H);
RX2 = reshape(rx,[winLen,1,nWins]);

tmp = pagefun(@mldivide,HH,Hp);
ws = mm3d(tmp,RX2); % weights

% figure(42)
% imagesc(abs(permute(ws,[3,1,2])))

direct = mm3d(H,ws); % the direct path 

rxClean = rx - direct(:);

rxClean([1:nTaps,end-nTaps:end]) = 0;


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