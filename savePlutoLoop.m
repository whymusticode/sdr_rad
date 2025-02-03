function savePlutoLoop(plutoID)

load params
sdr = sdrrx('Pluto', 'CenterFrequency', fc,'GainSource', 'Manual', ...
  'Gain', g,'BasebandSampleRate', fs,'OutputDataType',dataType,...
  'RadioID', ['usb:',num2str(plutoID)],'SamplesPerFrame',N);

ct = 0;
while 1
  tic
  ct = ct + 1;
  data = step(sdr);
  save(['.\data\',num2str(plutoID),'_',num2str(ct)])
  toc
end

