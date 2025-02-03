
% step 1: record data by running "a0" and "a1" simultaneously in 2
% seperate matlab instances
% step 2: align data, by loading both into one instance and then running
% "align.m"
% step 3: run this script to call "makeRDIs" to see RDIs

% a0
% a1
% align

p.fs = fs; % sample frequency
p.fc = fc; % center frequency
p.rMax = 2e5; % maximum range
p.vMax = 300; % maximum doppler
p.cpi = 3; % seconds to integrate

% pass = gpuArray(DATA(1:round(end/2),:));
% dat = gpuArray(DATA);
profile on
[r,v,RDIs] = makeRDIs(DATA,p);
profile viewer


figure(95)
filename = 'RDIs.gif';
for i = 1:size(RDIs,3)
  aRdi = abs(RDIs(:,:,i));
  hold off
  imagesc(r/1e3,v,20*log10(aRdi/mean(aRdi(:))))
  xlabel('Range Difference of Arrival (km)')
  ylabel('Doppler shift (m/s)')
  colorbar
  colormap(jet)
  caxis([10,30])
  font
  pause(.2)
  
  h = figure(95);
  im = frame2im(getframe(h));
  [imind,cm] = rgb2ind(im,256);
  if i == 1
    imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',5);
  else
    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2);
  end
end


% writerObj = VideoWriter('RDIs.mp4');
% writerObj.FrameRate = 6;
% open(writerObj);
%     writeVideo(writerObj, getframe(h));

% close(writerObj);


