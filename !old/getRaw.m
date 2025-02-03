function getRaw(sdr,nPerFrame,nFrames,jj)





raw = zeros(nPerFrame,nFrames,'int16');


for i = 1:nFrames
  gainHist(i) = sdr.TunerGain;
  raw(:,i) = step(sdr);
  sdr.TunerGain = gainRule(sdr.TunerGain,raw(:,i));
end
TIMESTOP = datenum(clock)*24*3600;
raw = int8(raw);

save(['./data/',num2str(jj)],'raw','TIMESTOP','gainHist')











