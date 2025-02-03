
% in PUTTY: serial, COM5, 115200
% .EnableBurstMode = true % this will allow you to save off more data

clear
fc = 104.5e6; % (hz)
fs =  .3e6; % (hz)
g = 30; % (db)
timePerStep = .5; % (s)
steps = 1e2;
nPluto = 4;
dataType = 'int16';

N = round(fs*timePerStep);
for i = 1:nPluto
  sdr{i} = sdrrx('Pluto', 'CenterFrequency', fc,'GainSource', 'Manual', ...
    'Gain', g,'BasebandSampleRate', fs,'OutputDataType',dataType,...
    'RadioID', ['usb:',num2str(i-1)],'SamplesPerFrame',N);
end
for j = 1:nPluto;step(sdr{j});end % start connection with all devices


data = zeros(N,steps,nPluto,dataType);
for i = 1:steps
  disp(i)
  tic
  for j = 1:nPluto
    data(:,i,j) = step(sdr{j});
  end
  toc
end


return

data = int16(data);
20*log10(max(abs([real(data(:,1,1));imag(data(:,1,1))]))/128)

tic 
load('test')
toc

data = data(:,1,1);

tic
save('test','data')
toc

return
% length(unique(real(tmp)))
rx = single(data(:,:,2));
tx = single(data(:,:,4));
[lag,TX,RX,MAX] = getLag(rx(:),tx(:),1);

dat = data;
f = fc + fs*freqBins(length(dat));
fftData = fftshift(20*log10(abs(fft(dat))));
figure(3)
plot(smpl(f,1e1)/1e6,smpl(fftData,1e1))
font



p.vMax = 600;
p.rMax = 1e7;
p.fs = fs;
p.fc = fc;
[range,velocity,image] = sledge(rx(:),tx(:),p);
figure(6)
imagesc(range,velocity,20*log10(abs(image)))


%% listen to the radio
% rx = data(:,:,1);
% tx = data(:,:,2);
% rx = rx(:);
% tx = tx(:);
% deFac = 96*25;
% nSteps = floor(length(rx)/deFac);
% nNew = nSteps*deFac;
% rx = reshape(rx(1:nNew),[deFac,nSteps]);
% tx = reshape(tx(1:nNew),[deFac,nSteps]);
% playRx(rx, fs)
% playRx(tx, fs)






% (fc + freqBins(200)*fs)/1e6
%% other notes



% tic
% step(sdr{1});
% toc
% tic
% step(sdr{2});
% toc
% 
% % takes about 9 seconds to set up new connection
% 
% 
% 1e7

