
bb = 4:.5:8;
ssnnrr = 10:5:40;
res = [];
val = [];
for i = 1:length(bb)
  for j = 1:length(ssnnrr)
    

B = 10^bb(i); % bandwidth hz
snr = 10^(ssnnrr(j)/20); % linear, in voltage

t0 = 1/B/7;


t = (-1:1)/B/2; % times the signal is sampled

fun = @(x,t)x(1)*(1-(pi*B*(t-x(2))).^2/6); % parabolic approx
x = [snr,t0];
CRB = getCrb(fun,x,t,ones(size(t)));
res(i,j) = sqrt(CRB(2,2))*3e8;
% plot(t*3e8,fun(x,t))
% grid on
% xlabel('range')
% title('sinc function approximation')
% disp(['range accuracy is: ',num2str(round(sqrt(CRB(2,2))*3e8)),' m'])
  end
  val(i) = interp1(res(i,:),ssnnrr,20);
end



hold off
imagesc(ssnnrr,bb,log10(res))
hold on
plot(val,bb,'k','LineWidth',6)
xlabel('SNR (db)')
ylabel('Log10(bandwidth)')
colorbar
legend('20 meters, size of a fighter jet')
title('Log10( RMS range accuracy (m) )')
set(gca,'fontSize',15)