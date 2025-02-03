

%% this looks at the difference between pulse compression and least squares

% clear

bins = 2e2; % number of range bins computed
len = 2e4; % number of samples integrated for
del = [20,30,40]; % true target location delays in samples
dp2sig = 1000; % direct path to target signal ratio
nsr = 0; % noise to signal ratio (linear units)
subDelay = .2; % extra delay 
Vel = 300;
fs = .3e6;
fc = 1e8;
%% make truth
c = 3e8;
CPI = len/fs;


wfm = randn(len,1) + 1i*randn(len,1); % random waveform
H = convmtx(wfm,bins); H = H(1:len,:);

del = del(del < bins);
truth = zeros(bins,1);
truth(del) = 1;

wfm2 = timeShift(wfm,subDelay); %  subsample delay (issue with hardware)
H2 = convmtx(wfm2,bins); H2 = H2(1:len,:);
noise = nsr*randn(size(H,1),1) + nsr*randn(size(H,1),1)*1i;
dopShift = exp(linspace(0,1,len)'*CPI*1i*2*pi/c*Vel*fc);
y = dp2sig*wfm2 + H2*truth.*dopShift + noise;
  

scale = 100 / max(abs([real(y(:));imag(y(:))]));
y = double(int8(y * scale))/scale; % simulate bit depth 

%% signal processing
% tic
% pc = H'*y/(wfm'*wfm); % pulse compression 
% toc
profile on 
rdi = LSsledge(wfm,y,50,15);
% ls = (H'*H)^-1*H'*y; % least squares (computationally expensive)
profile viewer

dpp = c/len*fs/fc; % doppler per pixel 
rpp = c/fs;

%% plotting 
figure(203)

imagesc((0:49)*rpp/1e3,(-15:15)*dpp,20*log10(abs(rdi')))
caxis([-10,20])
xlabel('range (km)')
ylabel('velocity (m/s)')

% subplot(2,1,1)
% plot(20*log10(abs(pc)))
% title('compression')
% xlabel('range bin')
% font
% subplot(2,1,2)
% plot(20*log10(abs(ls)))
% title('least squares')
% xlabel('range bin')
% font


% l = chol(H'*H);

