function saveDataLoop(sdrRx,saveName,PLOT,gLim)

i = 0;
while true
  i = i + 1;
  [DATA,~,lost] = step(sdrRx);
  if lost ~= 0;warning('lost samples in auxiliary');end
  check = length(unique(real(DATA)));
  disp(check)
  
  if check < 2000 % my version of automatic gain control
    sdrRx.Gain = sdrRx.Gain + 2;
  elseif check > 3500
    sdrRx.Gain = sdrRx.Gain - 2;
  end
  if sdrRx.Gain > gLim;sdrRx.Gain = gLim; end
  DATA = DATA - mean(DATA);
  
  save([saveName,num2str(i)],'DATA')
  
  if PLOT
  
  tp = 20*log10(faf(DATA));
  plot(smpl(tp,100))
  font
  drawnow
  end
  
end