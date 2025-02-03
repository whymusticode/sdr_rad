function [data,sdr,lost,check] = myStep(sdr)
[tmp1,~,lost] = step(sdr);
data = tmp1 - mean(tmp1);

check = length(unique(real(tmp1)));

disp(['unique bits: ',num2str(check)])
% gLim = 80;
% if check < 150 % my version of automatic gain control
%   sdr.TunerGain = sdr.TunerGain + 2;
% elseif check > 200
%   sdr.TunerGain = sdr.TunerGain - 2;
% end
% if sdr.TunerGain > gLim;sdr.TunerGain = gLim; end

