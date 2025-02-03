function [DATA,sdrRx,check] = loadDataOnce(sdrRx,PLOT,gLim)

% i = 0;
% while true
%   i = i + 1;
[DATA,~,lost] = step(sdrRx);
if lost ~= 0;warning('lost samples in auxiliary');end
check = length(unique(real(DATA)));
% disp(check)

if check < 150 % my version of automatic gain control
    sdrRx.TunerGain = sdrRx.TunerGain + 2;
elseif check > 200
    sdrRx.TunerGain = sdrRx.TunerGain - 2;
end
if sdrRx.TunerGain > gLim;sdrRx.TunerGain = gLim; end
DATA = DATA - mean(DATA);

if ~ PLOT;return;end

tp = 20*log10(faf(DATA));
plot(smpl(tp,100))
font
drawnow


% end