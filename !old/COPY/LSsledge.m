function rdi = LSsledge(tx,rx,nR,nD)
% this creates a least squares range doppler map 
% with nR range bins and +- nD doppler bins



n = length(tx);

TX = fft(tx);
H = [];
dVec = -nD:nD;
for i = dVec
  txr = ifft(circshift(TX,[i,0]));
  tmp = convmtx(txr,nR); tmp = tmp(1:n,:);
  H = [H,tmp];
end

H = gpuArray(single(H));
rx = gpuArray(single(rx));

ls = (H'*H)\(H')*rx; 
% ls = (H')*rx./(tx'*tx); % pulse compression


rxc = rx - H*ls;
disp(20*log10(abs(rxc'*rxc/(rx'*rx))))

nD2 = length(dVec);
rdi = reshape(ls,[nR,nD2]);