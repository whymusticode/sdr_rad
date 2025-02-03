function [ra,ve,rdi] = sledgeFast(rx,tx,p)
% needs windowing 
% p.fs = .2e6; % sample frequency
% p.fc = 1e8; % center frequency
% p.rMax = 1e5; % maximum range
% p.vMax = 600; % maximum doppler

% profile on



% create bins to be computed
c = 3e8;
n = length(rx); % number of samples, assumed length(rx) == length(tx)
cpi = n/p.fs; % integration time (s)
dMax = p.vMax/c*p.fc; % max frequency 
nRan = round(p.rMax/c*p.fs); % number of bins to go
DOPS = 0:1/cpi:dMax; % gives me frequency
DOPS = [-fliplr(DOPS),DOPS(2:end)];
nDop = length(DOPS);
ra = (1:nRan)*c/p.fs; % range bins
ve = DOPS*c/p.fc; % velocity bins


nB = 1000; % number of batches
Nb = n/nB; % length of batches, can't be less than number of doppler bins
% optimally this is an integer number of doppler bins and power of 2

% tx = gather(tx);
% rx = gather(rx);


P = (1:Nb)';
r = permute(1:nB,[1,3,2]);
l = 1:nRan;
ref = conj(tx);

rnb = r*Nb;
idx1 = rnb+P;
idx1(idx1 > n) = idx1(idx1 > n)-n;
sr = rx(idx1);

% this is the part that takes forever and will run us out of memory
idx2 = idx1 - l;
idx2(idx2 > n) = idx2(idx2 > n)-n;
idx2(idx2 < 1) = idx2(idx2 < 1)+n;
sref = ref(idx2);



ccfs = sum(sref.*sr,1);

% m is column vector 
rdi = permute(fft(ccfs,[],3),[3,2,1]);
rdi = [rdi(end-floor(nDop/2)+1:end,:);rdi(1:round(nDop/2),:)];

% m = DOPS'*cpi;
% rdi = sum(exp(-1i*2*pi*m.*rnb/n).*ccfs,3);


% figure(34)
% imagesc(ra,ve,log(abs(rdi)))


% db2 = (nDop-1)/2;

% n*nRan*nDop/1e9 % order pf calculations for regular convolution
% nRan*(n+n*log2(n))/1e9 % for direct fft 
% (2*n*log2(n) + nDop*(n+n*log2(n)))/1e9 % for correlation fft

% batches
% Nb*nB*nRan/1e9


% put tx and rx in frequency domain on GPU 
% if ~isa(rx,'gpuArray');rx = gpuArray(single(rx));end
% if ~isa(tx,'gpuArray');tx = gpuArray(single(tx));end
% TX = fft(conj(flipud(tx)));
% RX = fft(rx);
% 
% % This is for speed when memory isn't an issue
% TTX = zeros(n,nDop,'single','gpuArray');
% TTX(:,1) = circshift(TX,[-db2,0]);
% i = 1;
% while i < nDop
%   idx = i+1:2*i;
%   idx = idx(idx <= nDop);
%   TTX(:,idx) = circshift(TTX(:,idx-i),[i,0]);
%   i = i*2;
% end
% rdi = ifft(RX.*TTX);
% rdi = rdi(n/2+1:n/2+nRan,:)';

% 
% 
% % Limits operation size by gpuMemory
% maxi = floor(gpuMem*1e9/32/n); % 8 bytes in a complex signal
% rdi = zeros(dBins,rBins,'single','gpuArray');
% TTX = circshift(TX,[-db2,0]);
% rdiPre = ifft(RX.*TTX);
% rdi(1,:) = rdiPre(n-buf:n+rBins,:)';
% i = 1; % how much to circshift and length of window
% start = 1;
% while start + i <= dBins
%   TTXnew = circshift(TTX,[i,0]); % this takes forever 
%   rdiPre = ifft(RX.*TTXnew); 
%   rdi(start+1:start+i,:) = rdiPre(n-buf:n+rBins,:)';
%   
%   
%   start = start + i;
%   
%   if i < maxi
%     i = i*2;
%     TTX = [TTX,TTXnew];
%   end
%   
%   
%   rdi(1,:)
%   
% end
% rdi = rdi(n-buf:n+rBins,:)';



% slow memory limited processing
% rdi = zeros(dBins,rBins,'single','gpuArray');%,'single','gpuArray'
% rotIdx = -db2:db2;
% for i = 1:length(rotIdx)
%   TTX = circshift(TX,[rotIdx(i),0]);
%   rdiPre = ifft(RX.*TTX);
%   rdi(i,:) = rdiPre(n-buf+1:n+rBins)';
% end




% profile viewer
%% end of used code

% this memory operation takes up most of the time
% TTXp = convmtx(TX,dBins);
% db2 = (dBins-1)/2;
% TTX = TTXp(db2+1:end-db2,:);


% shiftindex = (1:dBins) -(dBins-1)/2 - 2;
% S = full(sparse(mod(shiftindex,2*n)+1,1:dBins,1,2*n,dBins));
% TTX = ifft(fft(TX).*fft(S),'symmetric');



% tmp = circshift(TX,[-(dBins-1)/2,0]);
% TTX = toeplitz(tmp,[tmp(1);fliplr(tmp(end-dBins+2:end))]);

% idx = (1:2*n)' + (1:dBins)-(dBins-1)/2-1;
% idx(idx < 1) = idx(idx < 1) + 2*n;
% idx(idx > 2*n) = idx(idx > 2*n) - 2*n;
% TTX = TX(idx);

% tic
% TTXp = arrayfun(@(k)circshift(TX,k),(1:dBins)-(dBins-1)/2-1,'uni',0);
% toc
% TTX = [TTXp{:}];

% TTX = zeros(2*n,dBins,'single','gpuArray');
% for i = 1:dBins
%   TTX(:,i) = circshift(TX,[i-(dBins-1)/2-1,0]);
% end

%% for windowing
% WIN = repmat(hamming(n),[1,dBins]);
% rep = rep.*WIN;

% DSHIFT = gpuArray(exp(1i*(1:n)'*DOPS*2*pi*cpi/n)); 

% large matrix gpuArray(...
% rep = bsxfun(@times,DSHIFT,tx);

% DSHIFT = gpuArray(exp(1i*(1:n)')*gpuArray(DOPS*2*pi*cpi/n)); 
% rdi = rdi(n+1:n+rBins,:)';
% r = (1:rBins)*c/2/p.fs;

% figure(1)
% aRDI = abs(rdi);
% hold off
% imagesc(r/1e3,v,20*log10(aRDI/median(aRDI(:))))
% xlabel('Range Difference of Arrival (km)')
% ylabel('Doppler shift (m/s)')
% colorbar
% % caxis([0,20])
% colormap(jet)


%% This shows how a time shift is created by a phase shift in freq
%  t=0:0.001:0.1-0.001;
%  Fs = 1e3;
%  freq1 = 100;
%  x1=sinc(2*pi*freq1*t);
%  Delay=1.2;
%  yp = fft(x1);
%  yp = yp(1:length(x1)/2+1);
%  f = 0:Fs/length(x1):500;
%  yp = yp.*exp(-1i*2*pi*f*Delay*(1/Fs));
%  yp = [yp conj(fliplr(yp(2:end-1)))];
%  rep = ifft(yp,'symmetric');
%  hold off
%  plot(t(1:100),sinc(2*pi*freq1*(t-Delay/Fs)),'b:');
%  hold on;
%  plot(t(1:100),rep(1:100),'r');

