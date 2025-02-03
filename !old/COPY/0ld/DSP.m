

% 2^20 (1 million) is smallest fft GPU performs faster than CPU
% gpu is only ~20 times faster when maxing out ram 


% reset(gpuDevice)
%% params
% fs = 160e6; % sample rate
% cpi = .5; % seconds integration
bw = .25e6; % bandwidth
p.fc = fc; % center frequency
p.rMax = 2e5; % maximum range


% v = 100.8;% target velocity
% r = 100e3; % target range
% SNR = 1e-2; % target SNR

% nCh = 7;
% cfs = rand(1,1,nCh)*.8-.4;
nCh = 1;
cfs = 0;

nTaps = 1e2; % number of taps in direct path canceler



c = 3e8;

%% preliminary stuff
ds = floor(fs/bw); % downsample rate


p.fs = fs/ds; % sample frequency
% n = cpi*fs;

data = [DATA(:,2,1);DATA(:,2,2)];

dConv = round(-cfs*n);

FILTER = fft(fir1(1e4,1/ds)',2*n);

% make bins
n0 = n/ds;
rpp = c/p.fs; % range per pixel
range = rpp:rpp:p.rMax; % range bins
nR = length(range); % batch length
nB = floor(n0/nR); % number of batches
n = nR*nB; % truncated integration
dpp = c/n*p.fs/p.fc; % doppler per pixel
velocity = (-floor(nB/2):ceil(nB/2)-1)*dpp; % velocity bins


pfa = 1e-10;
tX = 5; % train gates in X
tY = 5; % train gates in Y
gX = 3; % guard gates in X direction
gY = 3; % guard gates in Y direction
Nkern = (2*(tX+gX)+1) *(2*(tY+gY)+1) - (2*gY+1)*(2*gX+1);
kernel = ones((gY+tY)*2+1,(gX+tX)*2+1);
kernel(tY+1:tY+2*gY+1,tX+1:tX+2*gX+1) = 0;
kernel(tY+gY+1,tX+gX+1) = 1/(1-pfa^(-1/Nkern));
CFAR = fft2(kernel,nR,nB);

%% real time processing
profile on

%% downconvert + filter + downsample
DAT = fft(data);
dat = zeros(n0*2,1,nCh,'single');
for i = 1:nCh
  for j = 1:1
    if dConv(i)>0
      down = ifft([DAT(end-dConv(i)+1:end,j);DAT(1:end-dConv(i),j)].*FILTER);
    else
      down = ifft([DAT(-dConv(i)+1:end,j);DAT(1:-dConv(i),j)].*FILTER);
    end
    dat(:,j,i) = down(2:ds:end,:,:);
  end
end

% clear data down FILTER
p.rMax = 4e5;
[range,velocity,image] = sledgeBatch(dat,dat,p);
aimage = abs(image);
imagesc(range/1e3,velocity,20*log10(aimage))
xlim([0,100])
ylim([-900,900])
xlabel('range (km)')
ylabel('velocity (m/s)')
font
title('ambiguity function for FM signal')
max(aimage(:))/median(aimage(:))
%% direct path cancelation goes here 
% 2 is transmit, 1 is recieve

match = dat(nTaps/2:nTaps:end,1,:);
ref = reshape(dat(:,2,:),[n0/nTaps,nTaps,nCh]);
refp = permute(ref,[2,1,3]);
taps = mm3d(mm3d(pagefun(@inv,mm3d(refp,ref)),refp),match);

%% Integrate then CFAR:
% 
dat2 = reshape(dat(1:n,:,:),[nR,nB,2,nCh]);
SR = fft(dat2(:,:,1,:));
SREF = fft(flipud( conj(dat2(:,:,2,:)) )); %conj( 
corFun1d = ifft(SR.*(SREF + circshift(SREF,[0,1]))); % range 
image = fftshift(fft(corFun1d,[],2),2); % doppler 
image2 = image.*conj(image);
% imagesc(log(abs(image(:,:,1,1))))

% RDIfreq = conj(fftshift(SR.*(SREF + circshift(SREF,[0,1])),2));
% ^ conj(fft(x)/numel(x)) == ifft(x)
test = ifft2(fft2(image2).*CFAR);

% tmp = ifft2(fft2(image2(:,:,1,1).*CFAR));
% tmp(1:4,1) 
% test(1:4,1,1,1)

profile viewer
imagesc(real(test(:,:,1,1)) )



