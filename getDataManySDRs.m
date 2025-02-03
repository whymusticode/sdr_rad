

if exist('hardware','var')
%   clearvars -except hardware
else
  clear
  hardware = sdrinfo;
end

time = 3600; % seconds to record
nPerFrame = 3e5;
fs  = 3e5; %(hz) 225e3 < fs <= 300e3 or 900e3 < fs <= 3200e3
g = 30;

% FC = [98.5,107.9,89.7]*1e6;
% FC = [97.3,91.7,93.3]*1e6; % current best candidates

% freqs = 3e2./dataM2(:,3);

% for ii = 1:floor(length(freqs)/4)

% FC = [104.5,92.5,89.7,93.7]*1e6; % 99.5
% FC = [100.7,104.1,107.9,90.9]*1e6; % 99.5
FC = [103.3,106.7,88.9,93.3]*1e6;
% FC = round(freqs((1:4) + 4*(ii-1))*1e6);
%    1.0e+04 *
%
%    -4.5982    0.5428    0.0236
%     1.3686    3.4790    0.0152
%     1.2662   -2.7575    0.0170
%     2.3096    0.8091    0.0136
%     0.4123    2.1539    0.0173

nSDRs = length(hardware);
if length(FC) < ceil(nSDRs/2)
  error('2 channels per frequency')
end

nFrames = round(time*fs/ nPerFrame);

for i = 1:nSDRs
  sdrRx{i} = comm.SDRRTLReceiver(...
    'CenterFrequency', FC(ceil(i/2)), ...
    'SampleRate',      fs, ...
    'SamplesPerFrame', nPerFrame, ...
    'RadioAddress',   num2str(i-1),...
    'EnableTunerAGC',  false, ...
    'TunerGain',g);
  dat(i).x = zeros(nPerFrame,nFrames,'int16');
end

gainHist = zeros(nSDRs,nFrames);
timeStamp = gainHist;
count = zeros(nSDRs,1); % counts number of frames output from each SDR
TIMESTART = datenum(clock)*24*3600;
parfor i = 1:nSDRs
  for j = 1:nFrames
    tmp = step(sdrRx{i});
    dat(i).x(:,j) = tmp;
    timeStamp(i,j) = datenum(clock)*24*3600;
    gainHist(i,j) = sdrRx{i}.TunerGain;
%     sdrRx{i}.TunerGain = gainRule(sdrRx{i}.TunerGain,tmp);
  end
end
% for i = 1:nSDRs
%   release(sdrRx{i})
% end


% contains this much useable data: 
% min(timeStamp(:,end)) - max(timeStamp(:,1))



% 
% for i = 0:nSDRs/2-1
%   ch2 = single(dat(1+2*i).x(:));
%   ch1 = single(dat(2+2*i).x(:));
%   [rx,tx,~] = alignFun(ch1,ch2);
%   TX = gpuArray(tx(end-3e5+1:end));
%   RX = gpuArray(rx(end-3e5+1:end));
%   rxClean = removeDirectVect(200,30,TX,RX);
%   ca1 = 10*log10(abs(rxClean'*rxClean/(RX'*RX)));
%   rxClean = removeDirectVect(200,30,RX,TX);
%   ca2 = 10*log10(abs(rxClean'*rxClean/(RX'*RX)));
%   disp([min(ca1,ca2),FC(i+1)/1e6,i+1])
%   drawnow
% end

% end

return
% cm = min(count);
% for i = 1:nSDRs
%   tmp = single(dat(i).x(:,count(i)-cm+3:count(i)));
%   mn = mean(single(dat(i).x));
%   
%   dat(i).x = tmp(:) - mn;
% end

% disp(cancela)

return
% figure(1)
% plot(real(dat(i).x(:)))
% return
%
% tmp = real(dat(1).x);
% rows = 3e5;
% tMp = vec2mat(tmp,rows)';

% for i = 230:400
% figure(151)
% ttt = abs(real(dat(1).x((1:3e5)+3e5*i)));
% hist(ttt)
% pause(.03)
% end

shunck1 = max(abs(tMp));
shunck2 = median(abs(tMp));
figure(152)
x = 1:length(shunck1);
plot(x,shunck1,x,shunck2)

% save('3chanCollect','FC','dat','TIMESTAMP','-v7.3')




for i = 0:2
  [me(i+1).rx,me(i+1).tx,smallLag] = alignFun(dat(1+2*i).x,dat(2+2*i).x);
end


return

%% signal processing

tic
p.rMax = 2e5; % maximum range
p.vMax = 600; % maximum doppler
p.fs = 3e5;
for j = 1:100
  idx = (1:3e5) + 3e5*(j-1);
  for i = 1:3
    
    TX = gpuArray(me(i).tx(idx));
    RX =  gpuArray(me(i).rx(idx));
    rxClean = removeDirectVect(200,30,TX,RX);
    
    p.fc = FC(i);
    [res(i).r,res(i).v, res(i).rdi(:,:,j)] = sledgeBatch(TX,rxClean,p);
    
    cancela(i) = 10*log10(abs(rxClean'*rxClean/(RX'*RX)));
  end
  disp([j,  cancela])
  
  % chan(i+1).rx = rx;
  % chan(i+1).tx = tx;
end
toc % see how long imaging takes

for j = 1:100
  for i = 1:3
    subplot(1,3,i)
    im = abs(res(i).rdi(:,:,j));
    im = (im./median(im(:))).^2;
    [X,Y] = detect(im,1e-7);
    hold off
    x = res(i).r/1e3;
    y = res(i).v;
    imagesc(x,y,10*log10(im))
    hold on
    
    
    [x0,y0] = subP(x,y,X,Y,im);
    
    plot(x0,y0,'wo')
    caxis([12,30])
    
    det(j,i).x = [x0*1e3;y0];
    det(j,i).id = 1:length(x0);
    det(j,i).R = repmat(diag([1e3,4].^2),[1,1,length(x0)]);
  end
  pause(.01)
end

t = 1:100;
p.phi = 10; % m/s^3
p.model = 'ca_doppler';
p.nRemove = 4;
p.nFirmDeclare = 5;
p.useGpu = 1;
p.maxInit = 10;
tracks = assocFast(t,det(:,2),p);




for i = 1:numel(det)
  nDet(i) = size(det(i).x,2);
  
end






% x = 1:1e3;
% plot(x,real(RX(1:1e3)),x,real(rxClean(1:1e3)))

% plot(fftshift(abs(fft(dat(4).x))))


%       tmp = tmp-mean(tmp);
%       tmp = int8(tmp*100);


