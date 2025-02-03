function plotRoutine(rdi,r,v,RX,RXclean,TX)

figure(45)

subplot(1,2,2)
aRdi = abs(rdi);
hold off
imagesc(r/1e3,v,20*log10(aRdi/median(aRdi(:))))
xlabel('Range Difference of Arrival (km)')
ylabel('Doppler shift (m/s)')
% title(i)
%   colorbar
colormap(jet)
caxis([10,30])
font



FRX = fft(RX);
  FTX = fft(conj(flipud(TX)));
  FRXc = fft(RXclean);
% figure(46)
subplot(3,2,1)
plot(1:1e3,real(RX(1:1e3)),1:1e3,real(RXclean(1:1e3)))
title([num2str(round(10*log10(abs(RXclean'*RXclean) / ...
  abs(RX'*RX)))),' db'])
legend('RX','RX with direct cancelled')%,'Location','SouthOutside'

subplot(3,2,3)
d1freq = smpl(20*log10(fftshift(abs(FRX))),1000);
d2freq = smpl(20*log10(fftshift(abs(FTX))),1000);
d3freq = smpl(20*log10(fftshift(abs(FRXc))),1000);
plot(1:length(d1freq),d1freq,1:length(d2freq),...
  d2freq,1:length(d3freq),d3freq)
legend('RX spectrum','TX spectrum','clean RX')

subplot(3,2,5)
n = length(RX);
cpiLag = (1:n) - floor(n/2);
aval = abs(fftshift(ifft(FRX.*FTX)));
plot(smpl(cpiLag,1e3),smpl(aval,1e3))
title('alignment check')
drawnow