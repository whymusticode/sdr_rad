function SX = bootResampLoc(SX,hFunc,y,R)
yHat = feval(hFunc,SX);
W  = pdfCust(y-yHat,R);
ind = resampstr(W);
if sum(ind) == 0
  ind = 1:size(SX,2);
  disp('resample fail')
end
SX = SX(:,ind);
end