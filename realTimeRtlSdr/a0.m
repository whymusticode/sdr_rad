
% this stores tx data (radio address '00') (no resistor, slaved clock) 
% tends to be weaker than the other signal

% clear
load params
saveName = './data/rx';
saveDataLoop(sdrRx,saveName,0,gLim)

