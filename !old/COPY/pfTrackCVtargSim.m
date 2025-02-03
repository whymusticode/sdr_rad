% function pfTrackCVtargSim


% close all
clear
addpath(genpath('C:\Users\bestm\Desktop\MATLAB\FILTERING'))

ts = 2; % seconds integration time, also seconds between detections
nE = 1; % number of emitters
nT = 1; % number of targets in scene

phi = .2; % m/s^2 acceleration allowed of target
rErr = 1e2; % m    range accuracy
dErr = 5; % m/s   doppler accuracy
%% generate random positions + velocities
exyz = [(rand(nE,2)-.5)*1e5,7000+5000*rand(nE,1)]; % emitter locations
txyz = [(rand(nT,2)-.5)*1e5,5e3*ones(nT,1)]; % target
vxyz = [randn(nT,2)*100,zeros(nT,1)]; % target velocities
vInput = permute(vxyz,[3,2,1]);

R = diag([rErr,dErr].^2);

Q = diag([0,0,0,phi*ts*ones(1,2),0].^2);

hFunc = @(SX)permute(Hbistat(   permute(SX(1:3,:),[3,1,2]),...
  permute(SX(4:6,:),[3,1,2])  ,exyz   ),[2,3,1]);

nSamp = 10000;
SX = mvnrnd([txyz,vxyz]',diag([ones(1,3)*1e4,ones(1,3)*1e2]),nSamp)';
i = 0;
while abs(txyz+vxyz*ts*(i-1))<8e4
  i = i+1;
  curPos = permute(txyz+vxyz*ts*(i-1),[3,2,1]);
  y = permute(Hbistat(curPos,vInput, exyz),[2,1,3]); 
  
  
  SX = propSXlocal(SX,ts,Q);
  SX = bootResampLoc(SX,hFunc,y,R);
  
  plot(squeeze(curPos(:,1,:))/1e3,squeeze(curPos(:,2,:))/1e3,'+',...
    SX(1,:)/1e3,SX(2,:)/1e3,'.')
  font
  axis equal
  xlim([-1,1]*80)
  ylim([-1,1]*80)
  drawnow

  
end

% end
% 
% function SX = propSXlocal(SX,ts,Q)
% SX(1:3,:) = ts*SX(4:6,:);
% SX = SX + mvnrnd(zeros(size(SX,1),1),Q,size(SX,2))';
% end
% function SX = bootResampLoc(SX,hFunc,y,R)
% yHat = feval(hFunc,SX);
% W  = pdfCust(y-yHat,R);
% ind = resampstr(W);
% if sum(ind) == 0
% ind = 1:size(SX,2);
% disp('resample fail')
% end
% SX = SX(:,ind);
% end

% function [P,E] = pdfCust(DX,S)
% 
% %   if size(M,2) == 1
% %     DX = X-repmat(M,1,size(X,2));  
%     E = 0.5*sum(DX.*(S\DX),1);
%     d = size(S,1);
%     E = E + 0.5 * d * log(2*pi) + 0.5 * log(det(S));
%     P = exp(-E);
% %   elseif size(X,2) == 1
% %     DX = repmat(X,1,size(M,2))-M;  
% %     E = 0.5*sum(DX.*(S\DX),1);
% %     d = size(M,1);
% %     E = E + 0.5 * d * log(2*pi) + 0.5 * log(det(S));
% %     P = exp(-E);
% %   else
% %     DX = X-M;  
% %     E = 0.5*DX'*(S\DX);
% %     d = size(M,1);
% %     E = E + 0.5 * d * log(2*pi) + 0.5 * log(det(S));
% %     P = exp(-E);
% %   end
% end

% M = x_0;
% P = P_0;
% nSamp = 1000;
% SX = gauss_rnd(M,P,nSamp);
% 
% MM_BS = zeros(size(M,1),size(Y,2));
% PP_BS = zeros(size(M,1),size(M,1),size(Y,2));
% 
% % Filtering loop for bootstrap filter
% for k = 1:size(Y,2)
%    SX = ungm_f(SX,k) + gauss_rnd(0,1,size(SX,2));
%    W  = gauss_pdf(Y(:,k),ungm_h(SX,k),1);
%    ind = resampstr(W);
%    SX = SX(:,ind);
%    M = mean(SX);
%    P = var(SX);
%    MM_BS(:,k)   = M;
%    PP_BS(:,:,k) = P;    
% end



%% end of used shit

%     plot(squeeze(y(:,1,:))/1e3,squeeze(y(:,2,:)),'.')
%   font
%   xlim([0,2e2])
%   ylim([-900,900])


% create grid
% del = 1e3;
% delX = del;
% delY = del;
% delZ = del;
% x = -7e4:delX:7e4;
% y = x;
% z = 0:delZ:3e4;
% [X,Y,Z] = ndgrid(x,y,z);
% X = X(:);
% Y = Y(:);
% XYZ = [X,Y,Z(:)];
% R1 = sqrt(sum(permute(exyz,[3,2,1]).^2,2));
% R2 = sqrt(sum((permute(exyz,[3,2,1])-XYZ).^2,2));
% R3 = sqrt(sum(XYZ.^2,2));
% R = R2 + R3;
% % R = permute(R,[3,1,2]);
% R = gpuArray(single(permute(R,[3,1,2])));
% % tic
% % y = zeros(nE,2,nT,nSteps);

