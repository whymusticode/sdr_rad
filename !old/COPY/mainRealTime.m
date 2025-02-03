

% notes:
% 3e8/fc/4*39.4% bunny ears should be a little less than this size (inches)
% 3-10 km flight altitude common for airliners
reset(gpuDevice)
clear
addpath('data')
% set universal parameters for run
nTaps = 30; % number of taps in cancellation filter
winLen = 200; % length of cancellation window
p.rMax = 100e3; % maximum range
p.vMax = 600; % maximum doppler
p.cpi = 1; % seconds to integrate
fc = 97.3e6; % Center frequency (Hz) [88.5]
fs  = .3e6; %(hz) 225e3 < fs <= 300e3 or 900e3 < fs <= 3200e3
n = 375e3; % samples/frame  [maximum is 375e3]
g = 35; % amplification at the beginning (db?)
gLim = 80; % maximum requested amplification
buffer = 4*fs; % maximum number of samples channels might be off
makeGif = 0;
filename = [num2str(fc/1e5),'.gif'];
applyFilter = 0;
filterDS = 2; % downsample rate determines how much is filtered



%% end user parameters


sdrRx = comm.SDRRTLReceiver(...
  'CenterFrequency', fc, ...
  'EnableTunerAGC',  false, ...% false enables 'TunerGain'
  'SampleRate',      fs, ...
  'SamplesPerFrame', n, ...
  'TunerGain',g, ...
  'OutputDataType',  'single');

% save params for data loaders to use, then run the data loaders
% This opens up 2 new instances of MATLAB

% save params
% delete 'data\*.mat'
% !matlab -r a0 &
% !matlab -r a1 &

% return
nCPI = p.cpi*fs;
nCPI = round(nCPI/winLen)*winLen;% make CPI integer number of windows

% LAG = getLag(rx(1:length(tx)),tx,1);
% freqs = ((-floor(n/2):ceil(n/2)-1)/n*fs + fc)/1e6;

p.fs = fs;
p.fc = fc;

i = 1;
CPI = 0;
aux = 0; % approximate offset
LAG = 182209;
rx = [];
tx = [];
setLag = 0;
if applyFilter; filter = fft(fir1(1e3,1/filterDS)',nCPI);end
cpiLag = (1:nCPI) - floor(nCPI/2);
while true
  %   tic
  
  try % try to load the data
    load(['data\tx',num2str(i)])
    rx = [rx;DATA];
    load(['data\rx',num2str(i)])
    tx = [tx;DATA];
    
    if length(rx) ~= length(tx)
      disp('different lengths')
    end
    
  catch e % if the data isn't there yet, it likely doesn't exist
    disp(['pausing on step # ',num2str(i)])
    pause(n/fs)
    continue
  end
  i = i+1;
  
  if ~setLag
    
    if length(rx) <= buffer; continue; end
    % finding occasional length(rx) ~= length(tx)
    [val,lag] = xcorr(rx(end-buffer+1:end),tx(end-buffer+1:end));
    aval = abs(val);
    
    figure(145)
    plot(lag,aval)
    drawnow
    
    [VAL,loc] = max(aval);
    if ~(VAL/median(aval) > 50);continue;end
    LAG = lag(loc);
    disp(loc)
    setLag = true;
  end
  
  if length(rx) < abs(LAG) + nCPI; continue; end
  
  if LAG > 0;   D1 = rx(end-nCPI+1:end);  D2 = tx(end-nCPI+1-LAG:end-LAG);
  else;         D2 = tx(end-nCPI+1:end);  D1 = rx(end-nCPI+1+LAG:end+LAG);
  end
  
  rx = rx(end-nCPI+1:end);
  tx = tx(end-nCPI+1:end);
  
  if applyFilter
    D1 = ifft(fft(D1).*filter);
    D2 = ifft(fft(D2).*filter);
    
    D1 = D1(1:filterDS:end);
    D2 = D2(1:filterDS:end);
  end
  
  %   RX = D1; % to try processing on CPU  
  %   TX = D2;
  
  RX = gpuArray(D1);
  TX = gpuArray(D2);
  
  RX = RX - mean(RX);
  TX = TX - mean(TX);
  
  
  FRX = fft(RX);
  FTX = fft(conj(flipud(TX)));
  tmp = ifft(FRX.*FTX);
  aval = fftshift(abs(tmp));
  
  % make small corrections throughout 
  [~,loc] = max(aval);
  smallLag = cpiLag(loc);
  LAG = LAG + smallLag;
  if smallLag > 0; RX = RX(smallLag+1:end);  TX = TX(1:end-smallLag);
  else; RX = RX(1:end+smallLag);  TX = TX(-smallLag+1:end); end
  newLen = floor(length(TX)/winLen)*winLen;
  RX = RX(1:newLen);
  TX = TX(1:newLen);
  
  
  
%   tic

  RXclean = removeDirectVect(winLen,nTaps,TX,RX);
  RXclean = ifft(removeDirectVect(winLen,nTaps,fft(TX),fft(RXclean)));

  
  
  FRXc = fft(RXclean);
  
%   figure(403)
%   plot(fftshift(abs(FRXc)))
  
%   FRXc(abs(FRXc)>50) = 0;
%   RXclean = ifft(FRXc);
  
  
  RXclean([1:nTaps,end-nTaps:end]) = 0;
  [r,v, rdi] = sledgeBatch(TX,RXclean,p);
%   toc
  disp('CPI succesfully processed')
  CPI = CPI + 1;
  
  %% main plotting routine 
  figure(45)

  subplot(1,2,2)
  aRdi = abs(rdi);
  hold off
  imagesc(r/1e3,v,20*log10(aRdi/median(aRdi(:))))
  xlabel('Range Difference of Arrival (km)')
  ylabel('Doppler shift (m/s)')
  title(i)
  %   colorbar
  colormap(jet)
  caxis([12,30])
  
  
  subplot(3,2,1)
  plot(1:1e3,real(RX(1:1e3)),1:1e3,real(RXclean(1:1e3)))
  title([num2str(round(10*log10(abs(RXclean'*RXclean) / ...
    abs(RX'*RX)))),' db'])
  legend('RX','RX with direct cancelled')%,'Location','SouthOutside'
  
  subplot(3,2,3)
  d1freq = smpl(20*log10(fftshift(abs(FRX))),1000);
  d2freq = smpl(20*log10(fftshift(abs(FTX))),1000);
  d3freq = smpl(20*log10(fftshift(abs(FRXc))),1000);
  plot(1:length(d1freq),d1freq,1:length(d2freq),...
    d2freq,1:length(d3freq),d3freq)
  legend('RX spectrum','TX spectrum','clean RX')
  
  subplot(3,2,5)
  plot(cpiLag,aval)
  title('alignment check')
  drawnow
  
  if makeGif
    h = gcf;
    im = frame2im(getframe(h));
    [imind,cm] = rgb2ind(im,256);
    if CPI == 1
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',5);
    else
      imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2);
    end
  end
  
  
  
  %   toc
end



tic
getLag(rx,tx,1)
toc

%   RXclean(abs(RXclean)>.12) = 0; % a little clean up work



%   [DAT, ~,lost] = step(sdrRx);
%   if lost ~=0;warning('lost samples in auxiliary');end
%   check = length(unique(real(DAT)));
%   disp(check)
%   if check < 200 % my version of automatic gain control
%     sdrRx.TunerGain = sdrRx.TunerGain + 3;
%   elseif check > 250
%     sdrRx.TunerGain = sdrRx.TunerGain - 3;
%   end

% if lag is not 0, it is assumed everything is lined up

%
% D1 = gpuArray( ); % transmit
% D2 = gpuArray( ); % recieve
% rx = removeDirectVect(10,D2,D1);
% rx(abs(rx)>.12) = 0; % a little clean up work
% [r,v, rdi] = sledgeBatch(D1,rx,p);
%
% aRdi = abs(rdi);
% figure(95)
% hold off
% imagesc(r/1e3,v,20*log10(aRdi/mean(aRdi(:))))
% xlabel('Range Difference of Arrival (km)')
% ylabel('Doppler shift (m/s)')
% colorbar
% colormap(jet)
% caxis([10,30])
% font

%   tp = 20*log10(faf(DAT));
%   plot(smpl(tp,100))
%   font
%   drawnow
% end





% pass = gpuArray(DATA(1:round(end/2),:));
% dat = gpuArray(DATA);
% profile on
% [r,v,RDIs] = makeRDIs(DATA,p);
% profile viewer
%
% filename = 'RDIs.gif';
% % writerObj = VideoWriter('RDIs.mp4');
% % writerObj.FrameRate = 6;
% % open(writerObj);
% for i = 1:size(RDIs,3)
%
%
%   aRdi = abs(RDIs(:,:,i));
%
%   hold off
%   imagesc(r/1e3,v,20*log10(aRdi/mean(aRdi(:))))
%   xlabel('Range Difference of Arrival (km)')
%   ylabel('Doppler shift (m/s)')
%   colorbar
%   colormap(jet)
%   caxis([10,30])
%   font
% %   pause(.2)
%
%   h = figure(95);
%   im = frame2im(getframe(h));
%   [imind,cm] = rgb2ind(im,256);
%   if i == 1
%     imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',5);
%   else
%     imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2);
%   end
%
% %     writeVideo(writerObj, getframe(h));
%
% end
% close(writerObj);


