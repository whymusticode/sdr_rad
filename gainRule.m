function gOut = gainRule(gIn,tmp)



powMin = 60;
powMax = 80;
powGoal = 70;
gLim = 80;

tmp2 = abs(real(tmp));
powEst = median(tmp2);



if powEst > powMin && powEst < powMax
  gOut = gIn;
  return
end

change = 20*log10(powGoal/double(powEst));
% always shoot to get median(abs(real( == powGoal

gOut = gIn + change;

if gOut > gLim
  gOut = gLim;
end
% you will recieve the following error if output not a double:
% Changing the data type of a property value is not allowed without 
% first calling the release() method.





%% old way of doing it:
% check = length(unique(real(tmp)));
% 
% gLim = 60;
% 
% g = gIn;
% if check < 100 % my version of automatic gain control
%   g = g + 2;
% elseif check > 220
%   g = g - 2;
% end
% if g > gLim
%    g = gLim;
% end
% 
% 
% gOut = g;