function [r,v,RDIs] = makeRDIs(DATA,p)

%  writerObj = VideoWriter('RDIs.mp4');
%  writerObj.FrameRate = .2;
%   open(writerObj);
% ds = 2;
% FILTER = single(fft(fir1(1e4,1/ds)',n));

n = size(DATA,1); % length of data collected
nCPI = p.cpi*p.fs; % number of samples in a CPI 
i = 1;
while i*nCPI <= n
  D1 = gpuArray(DATA((i-1)*nCPI + 1:i*nCPI,1)); % transmit
  D2 = gpuArray(DATA((i-1)*nCPI + 1:i*nCPI,2)); % recieve
%   D1 = DATA((i-1)*nCPI + 1:i*nCPI,1);
%   D2 = DATA((i-1)*nCPI + 1:i*nCPI,2);
  
%   profile on
  rx = removeDirectVect(300,26,D2,D1);
%   profile viewer
%   return
% rx(abs(rx)>.12) = 0;
  
%     figure(95)
%     hold off
%     plot(real(rx(round(end/2):round(end/2)+1e3)))
%     hold on
%     plot(real(D2(round(end/2):round(end/2)+1e3)))
%     legend('canceled','direct')
    disp([num2str(round(10*log10(abs(rx'*rx) / abs(D2'*D2)))),' db'])
  
  [r,v, RDIs(:,:,i)] = sledgeBatch(D1,rx,p);
 
  
  i = i+1;
  
%   aRdi = abs(rdi);
%   figure(95)
%   hold off
%   imagesc(r/1e3,v,20*log10(aRdi/mean(aRdi(:))))
%   xlabel('Range Difference of Arrival (km)')
%   ylabel('Doppler shift (m/s)')
%   colorbar
%   colormap(jet)
%   caxis([10,30])
%   font
%   drawnow
%   
%   writeVideo(writerObj, getframe);
end


