
% This code proves out the ability to use >3 bistatic detection to solve
% ambiguity problems

% computation time scales with nT ^ nE
clear
addpath(genpath('C:\Users\bestm\Desktop\MATLAB\FILTERING'))

nSteps = 50;
nParticles = 1e4;
ts = .5; % seconds integration time, also seconds between detections
phi = 10; % m/s^3 process noise
nE = 7; % number of emitters
nT = 6; % number of targets in scene
rangeErr = 50; % [m] rms range error
dopErr = 10; % [m/s] rms doppler error
prob = 1; % probability of one detection getting through
reSeed = .05; % re - seed rate 
%% end user parameters
dims = 9;
FFF = [1,ts,ts^2/2;0,1,ts;0,0,1];
QQQ = [phi^2*ts^5/20,phi^2*ts^4/8,phi^2*ts^3/6;
  phi^2*ts^4/8,phi^2*ts^3/3,phi^2*ts^2/2;
  phi^2*ts^3/6,phi^2*ts^2/2,phi^2*ts];
z33 = zeros(3,3);
F = [FFF,z33,z33;z33,FFF,z33;z33,z33,FFF];
Q = [QQQ,z33,z33;z33,QQQ,z33;z33,z33,QQQ];

%% generate random positions
exyz = [(rand(nE,2)-.5)*1e5,300+100*randn(nE,1)]; % emitter locations
txyz = [(rand(nT,2)-.5)*1e5,5e3*ones(nT,1)+10e3*rand(nT,1)]; % target locations
% txyz = [(rand(nT,2)-.5)*1e3+2e4,10e3*ones(nT,1)]; % close targets
vxyz = [randn(nT,2)*200,zeros(nT,1)]; % target velocities
rex = [1,4,7,2,5,8,3,6,9];
X = zeros(dims,nSteps);
P = zeros(dims,dims,nSteps);
SX = [[(rand(nParticles,2)-.5)*1e5,5e3*ones(nParticles,1)+...
  10e3*rand(nParticles,1)]';...
  [randn(nParticles,2)*200,zeros(nParticles,1)]';zeros(3,nParticles)];
SX = SX(rex,:);
aCnt = ones(1,nParticles);

steps = 3;
redVal = 100*(100^(1/steps)).^-(0:steps-1);

tic
for j = 1:nSteps
  %% real world
  curPos = txyz+vxyz*ts*(j-1);
  y = Hbistat(permute(curPos,[3,2,1]),permute(vxyz,[3,2,1]),...
    exyz); % get all possible exact measurements
  
  idx = rand(nT,nE) < prob;
    
  nPerEmit = sum(idx,1);
  for i = 1:nE % put into realistic format with noise
    det(i).s = exyz(i,:);
    det(i).y = permute(y(i,:,idx(:,i)),[3,2,1]) + ...
      randn(nPerEmit(i),2).*[rangeErr,dopErr];
    det(i).R = repmat(diag([rangeErr,dopErr].^2),[1,1,nPerEmit(i)]);
  end
  
  %% filter 
  ind2 = rand(1,nParticles) < reSeed; % reseeding filter
  nres = sum(ind2);
  SX(rex,ind2) = [[(rand(nres,2)-.5)*1e5,5e3*ones(nres,1)+ 10e3*rand(nres,1)]';...
  [randn(nres,2)*200,zeros(nres,1)]';zeros(3,nres)];
  aCnt(ind2) = 1;

  Qnoise = gauss_rnd(zeros(dims,1),Q,size(SX,2));
  for i = 1:steps
    Qnoise(:,aCnt == i) = Qnoise(:,aCnt == i)*sqrt(redVal(i));
  end
  
  SX = F*SX + Qnoise; % process noise 
  
  OBS = Hbistat(permute(SX([1,4,7],:),[3,1,2]),...
    permute(SX([2,5,8],:),[3,1,2]),exyz);
  d2 = zeros(nE,nParticles);
  for i = 1:nE
    res = permute(det(i).y - OBS(i,:,:), [2,4,1,3]);
    resp = permute(res,[2,1,3,4]);
    si = inv2d(det(i).R);
    d2(i,:) = min(permute(mm3d(mm3d(resp,si),res),[3,4,1,2]));
  end
  d2Order = sort(d2);
  D2 = sum(d2Order(1:4,:));
  
  aCntTmp = aCnt;
  for i = 1:steps
    id = aCnt == i;
    W = chi2pdf(D2(id)/redVal(i),8);
    ind = resampstr(W);
    tmp = SX(:,id);
    SX(:,id) = tmp(:,ind);
    aCntTmp(id) = i+1;
  end
  
  id = aCnt > steps;  
  W = chi2pdf(D2(id),8);
  ind = resampstr(W);
  tmp = SX(:,id);
  aTmp = aCnt(id);
  try SX(:,id) = tmp(:,ind);aCntTmp(id) = aTmp(ind) + 1;catch; end
  
  aCnt = aCntTmp;
  
  
  
  figure(20)
  plot(curPos(:,1)/1e3,curPos(:,2)/1e3,'+',SX(1,:)/1e3,SX(4,:)/1e3,'.')
  xlim([-60,60])
  ylim([-60,60])
  font
  drawnow
  %   pause(.2)
end
toc

