% clearvars -except DATA fs fc n
% KERN = gpuArray(fir1(1e2,1/ds));
% DATA = gpuArray(DATA);
% data = conv(DATA(:,1),gather(KERN),'same');
% data2 = conv(DATA(:,2),gather(KERN),'same');
% DATA = gpuArray(DATA);

ds = 2;
FILTER = single(fft(fir1(1e4,1/ds)',n));


nCpi = size(DATA,3);
% n = size(DATA,1);
for i = 1:nCpi
  data = DATA(:,1,i);
  data2 = DATA(:,2,i);
  [val,lag] = xcorr(data,data2);
  figure(45)
  plot(lag,abs(val))
  [~,loc] = max(abs(val)); LAG = lag(loc);
  if LAG > 0;D1 = data(LAG:end); D2 = data2(1:end-LAG+1);
  else;LAG = -LAG; D2 = data2(LAG:end); D1 = data(1:end-LAG+1);
  end
%   [val,lag] = xcorr(D1,D2);
%   plot(lag,abs(val))


%   nNew = length(D1);
% FILTER = single(fft(fir1(1e4,1/ds)',nNew));
% D1 = ifft(fft(D1).*FILTER);
%   D2 = ifft(fft(D2).*FILTER);
%   D1 = D1(1:2:end);
%   D2 = D2(1:2:end);

%   tic
  rx = removeDirect(20,D2,D1);
%   figure(95)
%   hold off
%   plot(real(rx(round(end/2):round(end/2)+1e3)))
%   hold on
%   plot(real(D2(round(end/2):round(end/2)+1e3)))
%   legend('canceled','direct')
%   title([num2str(round(10*log10(abs(rx'*rx) / abs(D2'*D2)))),' db'])
%   toc
    
  
  
  p.fs = fs; % sample frequency
  p.fc = fc; % center frequency
  p.rMax = 1e5; % maximum range
  p.vMax = 200; % maximum doppler
%   rx = rxPre(1000:end);
%   tx = D1(1000:end);
%   tic



  [r,v,rdi] = sledgeBatch(D2,rx,p);
  
    aRdi = abs(rdi);
    figure(95)
    hold off
  imagesc(r/1e3,v,20*log10(aRdi/mean(aRdi(:))))
  xlabel('Range Difference of Arrival (km)')
  ylabel('Doppler shift (m/s)')
  colorbar
  colormap(jet)
  caxis([10,30])
  font
%   xlim([0,100])
  pause(.04)
%   store(:,:,i) = rdi;
%   toc
  % caxis([10,15])
  % drawnow
end
% return

%   plot(1:length(rx),real(tx),1:length(rx),real(rx))

% for i = 1:nCpi
%   aRdi = abs(store(:,:,i));
%   imagesc(r/1e3,v,20*log10(aRdi/mean(aRdi(:))))
%   xlabel('Range Difference of Arrival (km)')
%   ylabel('Doppler shift (m/s)')
%   colorbar
%   colormap(jet)
%   caxis([10,15])
%   xlim([0,100])
%   pause(.04)
% end

  % rxPre = D2 - opt';
  % plot(1:1e3,real(rxPre(1:1e3)),1:1e3,real(D2(1:1e3)))
  %
%   disp(10*log10(rxPre'*rxPre / (D2'*D2)))
  % rxPre = out(1e3:end);
  % n = length(D2);
  % plot(1:n,real(D2-LMSest'),1:n,real(D2-RLSest'),1:n,real(D2-opt'))
  %
  % RLSest(end-100:end)*RLSest(end-100:end)'
  % opt(end-100:end)*opt(end-100:end)'
  
%   plot(1:length(rx),real(rx),1:length(rx),real(tx))

% d1 = D1(1:ds:end);
% d2 = D2(1:ds:end);
% len = length(d1);
% 
% [val,lag] = xcorr(d1,d2);
% 
% plot(lag,abs(val))
% plot(1:len,real(d1),1:len,real(d2))


% h = randn(100,4)+1i*randn(100,4);
% x = randn(4,1) + 1i*randn(4,1);
% y = h*x;
%
%
% x2 = (h'*h)^-1*h'*y;



% o = 100; % (filter taps-1)/2
% 
% O = 2*o+1;

%% LS
% n = 1e5;
% H = zeros(n,O);
% y = zeros(n,1);
% for i = 1:n
%   y(i) = d2(i+o);
%   H(i,:) = d1(i:i+2*o)';
%   
% end

% y = real(y);
% H = real(H);
% X = (H'*H)^-1*H'*y;
% X = zeros(O,1);
% X(o+1) = 1;
% plot(1:n,real(y),1:n,real(H*X))
% plot(1:n,real(y),1:n,real(H(:,o+1)))
% res = y-H*X;
% disp(10*log10(res'*res/(y'*y)))
% return
%% RLS
% l = 1.01;
% I = eye(O);
% P = I*1e10;
% X = zeros(O,1);
% res = zeros(len,1);
% for i = o+1:1e3-o
%   y = d2(i);
%   H = d1(i-o:i+o)';
%
%
%   M = P*l;
%   S = H*M*H' + 1;
%   K = M*H'/S;
%   P = (I-K*H)*M;
%
%   res(i) = y-H*X;
%   X = X + K*(y-H*X);
% end
% plot(10*log10(abs(res)));
% font
% plot(abs(X))


%% RLS
% M = P*l;
% S = H*M*H' + 1;
% K = M*H'/S;
% P = (I-K*H)*M;
% X = X + K*(y-H*X);

%% kalman filter
% Xp = A*X;
% M = A'*P*A + Q;
% S = H*M*H' + R;
% K = M*H'*S^-1;
% P = (I-K*H)*M;
% X = Xp + K*(y - H*Xp);


