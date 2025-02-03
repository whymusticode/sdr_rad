
clear

cpi = .5; % seconds integration
p.fs = 8e6; % bandwidth
c = 3e8;
v = 100.8;%target velocity
p.fc = 600e6; % center frequency
r = 100e3; % range
SNR = 1e-1;
p.rMax = 2e5; % maximum range
% p.vMax = 600;
%% Create real data
n = round(p.fs*cpi); % number of samples in integration time
df = p.fc*v/c; % doppler HZ
sig = gpuArray(single(conv(randn(n,1) + ...
  1i*randn(n,1),fir1(100,.5),'same')));
noise = randn(n,1) + 1i*randn(n,1);
dShift = exp(1i*(1:n)'*2*pi*df*cpi/n);
rShift = [zeros(round(r/c*p.fs),1);sig(1:end-round(r/c*p.fs))];
RX = gpuArray(single(rShift.*dShift *SNR +noise)); % recieved signal


datPre = repmat([sig,RX],[1,1,7]);

%% Signal Processing
% profile on
tic
for i = 1:10
%   [RANGES,VELS,z2] = sledge(RX,sig,p);
  [RANGES,VELS,z2] = sledgeBatch(sig,RX,p);
end
toc
% profile viewer
% return

figure(1)
imagesc(RANGES/1e3,VELS,20*log10(abs(z2)/sqrt(n)/2))
xlabel('Range Difference of Arrival (km)')
ylabel('Doppler shift (m/s)')
colorbar
caxis([0,30])
colormap(jet)


%% end of code

return


% t = (1:rBins)/b;
% sig = exp(1i*pi*((1:n)'/n-.5).^2*n); % signal to pulse compress against
% plot(real(sig))


% dMax = vMax*2/c*fc;
% rBins = round(rMax*2/c*b); % number of bins to go
% DOPS = -dMax:1/(2*cpi):dMax;
% dBins =  length(DOPS);

% sig = conv(randn(n,1) + 1i*randn(n,1),fir1(100,.5),'same');

% RX = rShift.*dShift *SNR +noise; % recieved signal


% plot(abs(RX))

% disp(max(abs(z2(:)))/median(abs(z2(:))))

% DSHIFT = gpuArray(exp(1i*(1:n)'*DOPS*2*pi*cpi/n)); % large matrix
% rep = bsxfun(@times,DSHIFT,sig);
%
% x = repmat(RX,[1,dBins]);
% WIN = repmat(hamming(n),[1,dBins]);
% y = rep;%.*WIN;

% tic
% z2 = ifft(fft(x,2*n).*fft(conj(flipud(y)), 2*n));
% z2 = z2(n+1:n+rBins,:)';
% toc

% RANGES = (1:rBins)*c/2/b;
% VELS = DOPS*c/fc/2;

% a = 1:5;
% ifft(fft(a).*exp(1i*(0:4)/5))

%
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
%  y = ifft(yp,'symmetric');
%  hold off
%  plot(t(1:100),sinc(2*pi*freq1*(t-Delay/Fs)),'b:');
%  hold on;
%  plot(t(1:100),y(1:100),'r');

