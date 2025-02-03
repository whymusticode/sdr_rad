


% this code simulates bistatic range doppler detections from many emitters
% at a single reciever and attempts to disambiguate the detections to prove
% that the tracking algorithm functions as expected

% close all
clear

% addpath(genpath('C:\Users\bestm\Desktop\MATLAB\FILTERING'))

%% User parameters
PLOT = 1;
makeGif = 0;

plotTime = .1; % time between frames to watch plotting
ts = 2; % [s] integration time, time between detections
phi = 20; % [m/s^3] process noise (uncertainty in target motion)

nE = 7; % number of emitters

rangeErr = 1000; % [m] rms range error
dopErr = 20; % [m/s] rms doppler error

p.maxAcc = 50; % [m/s^2] maximum acceleration to initialize covariance
p.minR = 10e3; % [m] minimum detectable bistatic range
p.r00 = 1e5; % [m] r00 if the transmitters where monostatic radars
p.rcsMean = 0; % [dbsm] mean RCS
p.rcsStd = 5; % [dbsm] standard deviation of RCS
p.detectableSNR = 15; % SNR required for detection
p.lim = 5e4; % distance to end of the world

% pAssoc = .99; % probability of correct association

%% Tracking constant parameters
dims = 9;
FFF = [1,ts,ts^2/2;0,1,ts;0,0,1];
G = [ts^3/6,ts^2/2,ts];
QQQ = G'*phi^2*G;
F = zeros(9);
Q = F;
% thrBistat = chi2inv(pAssoc,2); % Mahalanobis distance criteria
thrBistat = 30;
for i = 0:2
  F([1,4,7]+i,[1,4,7]+i) = FFF;
  Q([1,4,7]+i,[1,4,7]+i) = QQQ;
end
p.F = F;
p.Q = Q;
p.thresh = thrBistat;
p.rErr = rangeErr;
p.dErr = dopErr;

%% make emitters
exyz = [(rand(nE,2)-.5)*2*p.lim,300+100*randn(nE,1)]; % emitter locations

%% Main loop
tStates = [genTarget(p),0,0,0]; % truth target states
det = struct([]); % detection data
tracks = det; % current tracks
TRACKS = det; % stored tracks that are dead
% profile on
iter = 0;
filename = 'myGif.gif';
meAn = [];
PP = [];
% delete(filename)
tic
while iter < 500
  iter = iter+1;
%   tic
  %% move targets and Simulate detection process
  [det,tStates,idx] = walkTruth(tStates,exyz,ts,p);
  %% tracker: propogation, association, update
  [det,tracks,TRACKS] = trackMaintain(det,tracks,TRACKS,p);
  
  %% track Initialization for unassociated detections
  [tracks,trNew] = trackInit(det,p,tracks);
  
  
  %% Plot routine
  if ~PLOT; continue;end
  figure(9)
  clf
  hold off
  detectable = sum(idx,2) > 3;
  plot(tStates(detectable,1)/1e3,tStates(detectable,2)/1e3,'g+')
  hold on
  plot(tStates(~detectable,1)/1e3,tStates(~detectable,2)/1e3,'r+')
  
  plotedTracks = false;
  for i = 1:length(tracks)
%     meAn = [meAn,tracks(i).xh(:,end)-tStates'];
%     PP = cat(3,PP,tracks(i).P(:,:,end));
    if size(tracks(i).xh,2) > 10
      plotedTracks = true;
    plot(tracks(i).xh(1,end)/1e3,tracks(i).xh(2,end)/1e3,'o','color',tracks(i).color)
%     quiver(tracks(i).xh(1,end)/1e3,tracks(i).xh(2,end)/1e3,...
%       tracks(i).xh(4,end)/3e1,tracks(i).xh(5,end)/3e1,'LineWidth',1)
    end
  end

  axis equal
  xlim([-1,1]*p.lim/1e3)
  ylim([-1,1]*p.lim/1e3)
  stepTime = toc;
  title(iter)
  font
%   pause(plotTime-toc)
  
  if makeGif
    frame = getframe(9);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if iter == 2
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
      imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
    end
  end
end
% profile viewer
toc
% plotTrackPerf(PP,meAn)


%% end of used code

%   legArg = {'targs with > 3 detections','other targets','tracks'};
%   legDx = [sum(detectable) ~=0 , sum(~detectable)~=0,plotedTracks];
%   legend(legArg(legDx),'Location','EastOutside')
  %   title(stepTime)


% for checking initialization
%   for i = 1:length(trNew)
%     res = tStates(:,1:6)'-trNew(i).xh(1:6);
%     m2 = sum(res.*(trNew(i).P(1:6,1:6)\res),1);
%     if ~any(m2<100)
%       disp('chicken')
%     end
%   end


% t = (0:nSteps-1)*ts;
%
%
%
% i = 1;
% err = [TRACKS(i).xh(1:3,:) - (txyz(i,1:3)'+...
%   t(1:size(TRACKS(i).xh,2)).*vxyz(i,:)') ; ...
%   vxyz(i,:)' - TRACKS(i).xh(4:6,:)];
% % plotTrackPerf(TRACKS(i).P,err,1:6)
% %
% figure(40)
% hold off
% plot(((txyz(:,1)+t.*vxyz(:,1))/1e3)',((txyz(:,2)+t.*vxyz(:,2))/1e3)')
% hold on
% for i = 1:length(TRACKS)
%   plot(TRACKS(i).xh(1,:)/1e3,TRACKS(i).xh(2,:)/1e3,':')
% end
% font



% %% generate random positions + velocities
% exyz = [(rand(nE,2)-.5)*1e5,300+100*randn(nE,1)]; % emitter locations
% txyz = [(rand(nT,2)-.5)*1e5,5e3*ones(nT,1)+10e3*rand(nT,1)]; % target
% % txyz = [(rand(nT,2)-.5)*1e3+2e4,10e3*ones(nT,1)]; % close targets
% vxyz = [randn(nT,2)*200,zeros(nT,1)]; % target velocities
% % rex = [1,4,7,2,5,8,3,6,9];
% X = zeros(dims,nSteps);
% P = zeros(dims,dims,nSteps);



% tracks = struct([]);
% for i = 1:nT
%   tracks(i).xh = [txyz(i,:),vxyz(i,:),zeros(1,3)]';
%   tracks(i).P = diag([ones(1,3)*rangeErr,...
%     ones(1,3)*dopErr,ones(1,3)*5].^2);
% end



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

