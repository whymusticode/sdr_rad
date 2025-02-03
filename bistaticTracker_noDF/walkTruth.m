
function [det,tStates,idx] = walkTruth(tStates,exyz,ts,p)

% currently this is SNR based, but it should be SIR based
% R1 = sqrt(sum(exyz.^2,2)); % interference dominates 

lim = p.lim; % extent on map
r00 = p.r00;
rcsStd = p.rcsStd;
rcsMean = p.rcsMean;
detectableSNR = p.detectableSNR;
rangeErr = p.rErr;
dopErr = p.dErr;

% walk states (constant velocity currently)
nT = size(tStates,1);

tStates(:,1:3) = tStates(:,1:3) + ts*tStates(:,4:6)+ ts^2/2*tStates(:,7:9);
tStates(:,4:6) = tStates(:,4:6) + ts*tStates(:,7:9);
tStates(:,7:8) = tStates(:,7:8) + 10*randn(nT,2);

tStates(:,4:6) = clip(tStates(:,4:6),400);
tStates(:,7:9) = clip(tStates(:,7:9),30);


% remove targets that have run off the screen
tStates = tStates( ~any(abs(tStates(:,1:2))>lim,2) ,:);

nE = size(exyz,1);

% randomly generate new targets
if nT < 15 && rand()<.2 % make 1 target
  x = genTarget(p);
  tStates = [tStates;x,0,0,0];
end

nT = size(tStates,1);

pos = permute(tStates(:,1:3),[3,2,1]);
vel = permute(tStates(:,4:6),[3,2,1]);

y = Hbistat(pos,vel, exyz); % possible measurements

R2 = sqrt(sum(pos.^2,2));
R3 = sqrt(sum((pos-exyz).^2,2));

idx = permute(10*log10(r00^4./(R2.*R3).^2),[3,1,2]) + ...
  randn(nT,nE)*rcsStd + rcsMean > detectableSNR;

R1 = sqrt(sum(exyz.^2,2));
% minRreq = permute(R2+R3-R1,[3,1,2]) > p.minR;
% disp(sum(idx(:))/numel(idx))

% idx = true(size(idx)); % just let all detections through

nPerEmit = sum(idx,1);
for j = 1:nE % put into realistic format with noise
  det(j).s = exyz(j,:);
  det(j).y = permute(y(j,:,idx(:,j)),[3,2,1]) + ...
    randn(nPerEmit(j),2).*[rangeErr,dopErr];
  det(j).R = repmat(diag([rangeErr,dopErr].^2),[1,1,nPerEmit(j)]);
end










