

% FC = 87.9:.2:107.9;
FC = 91.9;
p.fs = 250e3;
p.vMax = 600;
p.rMax = 200e3;
nCpi = 5e5;
nCh = numel(FC);

% data = store;
for j = 1:500
  rx = gS(data((1:nCpi)+nCpi*(j-1),2));
  tx = gS(data((1:nCpi)+nCpi*(j-1),1));
  rx0 = rx'*rx;
  
  rx = removeDirectVect(100,10,tx,rx);
  rx = ifft(removeDirectVect(100,5,fft(tx),fft(rx)));
%   rx = removeDirectVect(100,5,tx,rx);
  
  perf = 10*log10(real(rx'*rx./rx0)); % cancellation performance

  disp(perf)
%   if perf(i)> -10; continue; end
  
  p.fc = FC*1e6;
  [range,velocity,image] = sledgeBatch(tx,rx,p);
  imdb = 10*log10(real(image.*conj(image)));
  imdb=  imdb - median(imdb(:));
  imagesc(range/1e3,velocity,imdb,[10,30])
  title(FC)
  font
  drawnow
%   pause(2)
end



% disp([FC;perf]')






