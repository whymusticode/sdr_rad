
close all
clear
addpath(genpath('C:\Users\bestm\Desktop\MATLAB\FILTERING'))

nSteps = 100; % number of steps taken in simulation
ts = .5; % seconds integration time, also seconds between detections
phi = 10; % m/s^3 process noise
nE = 7; % number of emitters
nT = 6; % number of targets in scene
rangeErr = 50; % [m] rms range error
dopErr = 10; % [m/s] rms doppler error
reSeed = .05; % re - seed rate
r00 = 12e4; % r00 if the transmitters where monostatic radars
RCSstd = 5; % in db standard deviation of RCS
maxAcc = 50; % m/s^2 maximum acceleration to initialize covariance

thrBistat = chi2inv(.99,2);

%% end user parameters
dims = 9;
FFF = [1,ts,ts^2/2;0,1,ts;0,0,1];
QQQ = [phi^2*ts^5/20,phi^2*ts^4/8,phi^2*ts^3/6;
  phi^2*ts^4/8,phi^2*ts^3/3,phi^2*ts^2/2;
  phi^2*ts^3/6,phi^2*ts^2/2,phi^2*ts];
F = zeros(9);
Q = F;

for i = 0:2
  F([1,4,7]+i,[1,4,7]+i) = FFF;
  Q([1,4,7]+i,[1,4,7]+i) = QQQ;
end

%% generate random positions + velocities
exyz = [(rand(nE,2)-.5)*1e5,300+100*randn(nE,1)]; % emitter locations
txyz = [(rand(nT,2)-.5)*1e5,5e3*ones(nT,1)+10e3*rand(nT,1)]; % target
% txyz = [(rand(nT,2)-.5)*1e3+2e4,10e3*ones(nT,1)]; % close targets
vxyz = [randn(nT,2)*200,zeros(nT,1)]; % target velocities
% rex = [1,4,7,2,5,8,3,6,9];
X = zeros(dims,nSteps);
P = zeros(dims,dims,nSteps);

R1 = sqrt(sum(exyz.^2,2));
vInput = permute(vxyz,[3,2,1]);


% tracks = struct([]);
for i = 1:nT
  tracks(i).xh = [txyz(i,:),vxyz(i,:),zeros(1,3)]';
  tracks(i).P = diag([ones(1,3)*rangeErr*2,...
    ones(1,3)*dopErr*2,ones(1,3)*5].^2);
end

tic
det = struct([]);
TRACKS = det;
for i = 2:nSteps
  %% real world
  curPos = permute(txyz+vxyz*ts*(i-1),[3,2,1]);
  y = Hbistat(curPos,vInput, exyz); % get possible measurements
  
  R2 = sqrt(sum(curPos.^2,2));
  R3 = sqrt(sum((curPos-exyz).^2,2));
  
  idx = permute(10*log10(r00^4./(R2.*R3).^2),[3,1,2]) + ...
    randn(nT,nE)*RCSstd > 15; % 15 db SNR for detection
  
  nPerEmit = sum(idx,1);
  for j = 1:nE % put into realistic format with noise
    det(j).s = exyz(j,:);
    det(j).y = permute(y(j,:,idx(:,j)),[3,2,1]) + ...
      randn(nPerEmit(j),2).*[rangeErr,dopErr];
    det(j).R = repmat(diag([rangeErr,dopErr].^2),[1,1,nPerEmit(j)]);
  end
  
  
  
  %% tracker
  remove = false(size(tracks));
  for j = 1:length(tracks)
    xp = F*tracks(j).xh(:,end);
    M = F*tracks(j).P(:,:,end)*F' + Q;
    [yEst,H] = Hjac(xp,exyz);
    
    takenX = [];
    takenP = [];
    %     SS = [];
    HH = [];
    takeDx = false(nE*2,1);
%     count = 0;
    for k = 1:nE
      hh = [H([k,k+nE],:),zeros(2,3)];
      yy = yEst([k,k+nE]);
      S = det(k).R + hh*M*hh'; % residual covariance
      res = permute(det(k).y,[2,3,1]) - yy;
      M2 = mm3d(mm3d(permute(res,[2,1,3]),inv2d(S)),res);
      
      take = squeeze(M2) < thrBistat;
      if sum(take) > 1
%         return
%         disp('combined multiple');
      [~,take] = min(squeeze(M2));
      
      end
      if ~any(take); continue;end
%       count = count + 1;
      K = M*hh'*S(:,:,take)^-1;
      xp = xp + K*res(:,:,take);
      M = (eye(9) - K*hh)*M;
    end
    tracks(j).xh(:,end+1) = xp;
    tracks(j).P(:,:,end+1) = M;    
    
    if trace(M) > 1e3^2
      if isempty(TRACKS)
        TRACKS = tracks(j);
      else
      TRACKS(end+1) = tracks(j);
      end
      remove(j) = true;
    end
    
  end % loop over tracks
  tracks(remove) = [];
  
  
  
end
if isempty(TRACKS)
  TRACKS = tracks;
else
TRACKS(end+1:end+length(tracks)) = tracks;
end
toc
t = (0:nSteps-1)*ts;



i = 1;
err = [TRACKS(i).xh(1:3,:) - (txyz(i,1:3)'+...
  t(1:size(TRACKS(i).xh,2)).*vxyz(i,:)') ; ...
  vxyz(i,:)' - TRACKS(i).xh(4:6,:)];
% plotTrackPerf(TRACKS(i).P,err,1:6)
% 
figure(1)
hold off
plot(((txyz(:,1)+t.*vxyz(:,1))/1e3)',((txyz(:,2)+t.*vxyz(:,2))/1e3)')
hold on
for i = 1%:length(TRACKS)
  plot(TRACKS(i).xh(1,:)/1e3,TRACKS(i).xh(2,:)/1e3,':')
end
font


    %     x = iterSolve(S,r,d,rSig,dSig,x0)

    %       takeDx([k,k+nE]) = true;
    %       takenX = [takenX,det(k).y(take,:)'];
    %       takenP(end+1:end+2,end+1:end+2) = S(:,:,take);
    %       HH = [HH;hh];
    %     end
    %     HHH = [HH,zeros(size(HH,1),3)];
    %     HHH = [HHH(1:2:end,:);HHH(2:2:end,:)];
    %     K = M*HHH'*takenP^-1;
    %     fgfg = takenX';
    %     RES = fgfg(:) - yEst(takeDx);

%         takenP = [takenP,sqrt(diag(det(k).R(:,:,take)))];
%         SS = [SS;det(k).s];

