function [det,tracks,TRACKS] = assStep(det,tracks,TRACKS,trPar)
% hard association for bistatic range and doppler measurement only

F = trPar.F;
Q = trPar.Q;
thrBistat = trPar.thresh;

% track kill condition
tKC = [ones(3,1)*2e3;ones(3,1)*1e2;ones(3,1)*2e3];

nE = length(det);
exyz = zeros(nE,3);
for i = 1:nE
  exyz(i,:) = det(i).s;
  
  det(i).used = false(size(det(i).R,3),1);
end

remove = false(size(tracks));
for j = 1:length(tracks) % for every track 
  xp = F*tracks(j).xh(:,end); % update mean
  M = F*tracks(j).P(:,:,end)*F' + Q; % update covariance 
  [yEst,Hpre] = Hjac(xp,exyz); % determine "would be" observation 
  
  for k = 1:nE % consider every emitter
    H = [Hpre([k,k+nE],:),zeros(2,3)];
    yp = yEst([k,k+nE]); % y projected 
    S = det(k).R + H*M*H'; % residual covariance
    res = permute(det(k).y,[2,3,1]) - yp;
    M2 = mm3d(mm3d(permute(res,[2,1,3]),inv2d(S)),res);
    
    take = squeeze(M2) < thrBistat;
    
    % the case of multiple associated detections warants more consideration
    % in our application. 
    if sum(take) > 1  
      [~,take] = min(squeeze(M2));
    end
    
    if ~any(take); continue;end
    
    det(k).used(take) = true;
    
    K = M*H'*S(:,:,take)^-1;
    xp = xp + K*res(:,:,take);
    M = (eye(9) - K*H)*M;
  end
  
  tracks(j).xh(:,end+1) = xp;
  tracks(j).P(:,:,end+1) = M;
  
  
  % condition to kill a track
  
  tmp = sqrt(diag(M));
  if sum(tmp(4:5)) > 100
    if isempty(TRACKS)
      TRACKS = tracks(j);
    else
      TRACKS(end+1) = tracks(j);
    end
    remove(j) = true;
  end
  
end % loop over tracks
tracks(remove) = [];

% return remaining detections
for i = 1:nE  
  det(i).R = det(i).R(:,:,~det(i).used);
  det(i).y = det(i).y(~det(i).used,:);
end






