
% close all
clear
addpath(genpath('C:\Users\bestm\Desktop\MATLAB\FILTERING'))

nSteps = 100; % number of steps taken in simulation
ts = .5; % seconds integration time, also seconds between detections
nE = 3; % number of emitters
nT = 5; % number of targets in scene
% phi = 10; % m/s^3 process noise
% rangeErr = 50; % [m] rms range error
% dopErr = 10; % [m/s] rms doppler error
% r00 = 12e4; % r00 if the transmitters where monostatic radars
% RCSstd = 5; % in db standard deviation of RCS
% maxAcc = 50; % m/s^2 maximum acceleration to initialize covariance
% snrDetect = 15; % db for detect

%% generate random positions + velocities
exyz = [(rand(nE,2)-.5)*1e5,7000+5000*rand(nE,1)]; % emitter locations
txyz = [(rand(nT,2)-.5)*1e5,5e3*ones(nT,1)+10e3*rand(nT,1)]; % target
% txyz = [(rand(nT,2)-.5)*1e3+2e4,10e3*ones(nT,1)]; % close targets
vxyz = [randn(nT,2)*100,zeros(nT,1)]; % target velocities


% R1 = sqrt(sum(exyz.^2,2));
vInput = permute(vxyz,[3,2,1]);

% det = struct([]);

% create grid
del = 1e3;
delX = del;
delY = del;
delZ = del;
x = -7e4:delX:7e4;
y = x;
z = 0:delZ:3e4;
[X,Y,Z] = ndgrid(x,y,z);
X = X(:);
Y = Y(:);
XYZ = [X,Y,Z(:)];
R1 = sqrt(sum(permute(exyz,[3,2,1]).^2,2));
R2 = sqrt(sum((permute(exyz,[3,2,1])-XYZ).^2,2));
R3 = sqrt(sum(XYZ.^2,2));
R = R2 + R3;
% R = permute(R,[3,1,2]);
R = gpuArray(single(permute(R,[3,1,2])));
nBins = size(R,2);
% tic
% y = zeros(nE,2,nT,nSteps);
% err = zeros(nE,nBins);
for i = 1:nSteps
  %% real world
  curPos = permute(txyz+vxyz*ts*(i-1),[3,2,1]);
  y = Hbistat(curPos,vInput, exyz); %(:,:,:,i)
  
% err(:) = 0;
  err = abs(y(:,1,:) - R) <= del/2;
  
  num = sum(sum(err,1),3) == nE;
  
  plot3(X(num)/1e3,Y(num)/1e3,Z(num)/1e3,'.',...
    squeeze(curPos(:,1,:))/1e3,squeeze(curPos(:,2,:))/1e3,...
    squeeze(curPos(:,3,:))/1e3,'+')
    font
    axis equal
    axis([[-1,1]*80,[-1,1]*80,0,30])
  
%   
  drawnow
%     pause(.05)

  
end
% toc


%     plot(squeeze(y(:,1,:))/1e3,squeeze(y(:,2,:)),'.')
%   font
%   xlim([0,2e2])
%   ylim([-900,900])
 
  % ^ get possible measurements
  
%   R2 = sqrt(sum(curPos.^2,2));
%   R3 = sqrt(sum((curPos-exyz).^2,2));
%   
%   idx = permute(10*log10(r00^4./(R2.*R3).^2),[3,1,2]) + ...
%     randn(nT,nE)*RCSstd > snrDetect; % 15 db SNR for detection
  

%   nPerEmit = sum(idx,1);
%   for j = 1:nE % put into realistic format with noise
%     det(j).s = exyz(j,:);
%     det(j).y = permute(y(j,:,idx(:,j)),[3,2,1]) + ...
%       randn(nPerEmit(j),2).*[rangeErr,dopErr];
%     det(j).R = repmat(diag([rangeErr,dopErr].^2),[1,1,nPerEmit(j)]);
%   end
  
%   plot(squeeze(y(:,1,:,i))/1e3,squeeze(y(:,2,:,i)),'.')
%   font
%   xlim([0,2e2])
%   ylim([-900,900])
  

% t = (1:nSteps)*ts;
% plot(t,squeeze(y(1,2,:,:)))




% rex = [1,4,7,2,5,8,3,6,9];
% X = zeros(dims,nSteps);
% P = zeros(dims,dims,nSteps);
% thrBistat = chi2inv(.99,2);

%% end user parameters
% dims = 9;
% FFF = [1,ts,ts^2/2;0,1,ts;0,0,1];
% QQQ = [phi^2*ts^5/20,phi^2*ts^4/8,phi^2*ts^3/6;
%   phi^2*ts^4/8,phi^2*ts^3/3,phi^2*ts^2/2;
%   phi^2*ts^3/6,phi^2*ts^2/2,phi^2*ts];
% F = zeros(9);
% Q = F;
% 
% for i = 0:2
%   F([1,4,7]+i,[1,4,7]+i) = FFF;
%   Q([1,4,7]+i,[1,4,7]+i) = QQQ;
% end


% tracks = struct([]);
% for i = 1:nT
%   tracks(i).xh = [txyz(i,:),vxyz(i,:),zeros(1,3)]';
%   tracks(i).P = diag([ones(1,3)*rangeErr*2,...
%     ones(1,3)*dopErr*2,ones(1,3)*5].^2);
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

