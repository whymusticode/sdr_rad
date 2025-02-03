nCPI = p.cpi*fs;
nCPI = round(nCPI/winLen)*winLen;% make CPI integer number of windows

p.fs = fs;
p.fc = fc;

tuneF = [98.5,99.5,100.7,101.3];
fcs = (tuneF*1e6 - fc)/fs;

i = 1;
CPI = 0;
aux = 0; % approximate offset
LAG = -746497;
rx = [];
tx = [];
setLag = 1;
if applyFilter; filter = fft(fir1(1e3,1/filterDS)',nCPI);end
% while true
%   %   tic
%   
% %   try % try to load the data
%     load(['data\rx',num2str(i)])
%     rx = [rx;DATA];
%     load(['data\tx',num2str(i)])
%     tx = [tx;DATA];
% %   catch e % if the data isn't there yet, it likely doesn't exist
% %     disp(['pausing on step # ',num2str(i)])
% %     pause(n/fs)
% %     continue
% %   end
%   i = i+1;
%   disp(i)
% end

for i = 1:10
  cpiLag = (1:nCPI) - floor(nCPI/2);
  lag = getLag(rxtx((1:nCPI) + i*nCPI,1),rxtx((1:nCPI),2));
  drawnow
end

% return
  
%   
%   if ~setLag
%     
% %     if length(rx) <= buffer; continue; end
%     % finding occasional length(rx) ~= length(tx)
%     [val,lag] = xcorr(rx(end-buffer+1:end),tx(end-buffer+1:end));
%     aval = abs(val);
%     
%     subplot(2,2,1)
%     plot(lag,aval)
%     drawnow
%     
%     [VAL,loc] = max(aval);
%     if ~(VAL/median(aval) > 50);continue;end
%     LAG = lag(loc);
%     disp(loc)
%     setLag = true;
%   end
%   
%   if length(rx) < abs(LAG) + nCPI; continue; end
%   
%   if LAG > 0;   D1 = rx(end-nCPI+1:end);  D2 = tx(end-nCPI+1-LAG:end-LAG);
%   else;         D2 = tx(end-nCPI+1:end);  D1 = rx(end-nCPI+1+LAG:end+LAG);
%   end
%   
%   rx = rx(end-nCPI+1:end);
%   tx = tx(end-nCPI+1:end);
%   
%   if applyFilter
%     D1 = ifft(fft(D1).*filter);
%     D2 = ifft(fft(D2).*filter);
%     
%     D1 = D1(1:filterDS:end);
%     D2 = D2(1:filterDS:end);
%   end
%   
%   %   RX = D1; % to try processing on CPU (bad idea) 
%   %   TX = D2;
%   
%   RX = gpuArray(D1);
%   TX = gpuArray(D2);
%   
%   RX = RX - mean(RX);
%   TX = TX - mean(TX);
%   
%   
%   FRX = fft(RX);
%   FTX = fft(conj(flipud(TX)));
%   tmp = ifft(FRX.*FTX);
%   aval = fftshift(abs(tmp));
% %   [~,loc] = max(aval);
% %   smallLag = cpiLag(loc);
% %   LAG = LAG + smallLag;
%   
% %   if smallLag > 0; RX = RX(smallLag+1:end);  TX = TX(1:end-smallLag);
% %   else; RX = RX(1:end+smallLag);  TX = TX(-smallLag+1:end); end
% %   
% %   newLen = floor(length(TX)/winLen)*winLen;
% %   RX = RX(1:newLen);
% %   TX = TX(1:newLen);
% %   tic
%   RXclean = removeDirectVect(winLen,nTaps,TX,RX);
%   
%   RXclean([1:nTaps,end-nTaps:end]) = 0;
%   [r,v, rdi] = sledgeBatch(TX,RXclean,p);
% %   toc
%   disp('CPI succesfully processed')
%   CPI = CPI + 1;
%   
%   %% main plotting routine 
%   figure(45)
% 
%   subplot(1,2,2)
%   aRdi = abs(rdi);
%   hold off
%   imagesc(r/1e3,v,20*log10(aRdi/median(aRdi(:))))
%   xlabel('Range Difference of Arrival (km)')
%   ylabel('Doppler shift (m/s)')
%   title(i)
%   %   colorbar
%   colormap(jet)
%   caxis([12,30])
%   
%   
%   subplot(3,2,1)
%   plot(1:1e3,real(RX(1:1e3)),1:1e3,real(RXclean(1:1e3)))
%   title([num2str(round(10*log10(abs(RXclean'*RXclean) / ...
%     abs(RX'*RX)))),' db'])
%   legend('RX','RX with direct cancelled')%,'Location','SouthOutside'
%   
%   subplot(3,2,3)
%   d1freq = smpl(20*log10(fftshift(abs(FRX))),1000);
%   d2freq = smpl(20*log10(fftshift(abs(FTX))),1000);
%   plot(1:length(d1freq),d1freq,1:length(d2freq),d2freq)
%   legend('RX spectrum','TX spectrum')
%   
%   subplot(3,2,5)
%   plot(cpiLag,aval)
%   title('alignment check')
%   drawnow
%   
%   if makeGif
%     h = gcf;
%     im = frame2im(getframe(h));
%     [imind,cm] = rgb2ind(im,256);
%     if CPI == 1
%       imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',5);
%     else
%       imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2);
%     end
%   end
%   
%   
%   
%   %   toc
% % end
% 
% 
% 
% 
% 
% 
% 
% 
% 
