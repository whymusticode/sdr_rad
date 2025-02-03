function [det,tracks,TRACKS] = trackMaintain(det,tracks,TRACKS,trPar)
% hard association for bistatic range and doppler measurement only

nTr = numel(tracks);
if nTr == 0
  return
end

rErr = trPar.rErr;
dErr = trPar.dErr;

F = trPar.F; % mean propagation
Q = trPar.Q; % covariance expansion
thrBistat = trPar.thresh; % m2 threshold

nE = length(det);
exyz = zeros(nE,3);
for i = 1:nE
  exyz(i,:) = det(i).s;
  det(i).used = false(size(det(i).R,3),1);
end

% update and put all the tracks into bistatic space for each emitter
for i = 1:length(tracks) % for every track 
  tracks(i).xh(:,end+1) = F*tracks(i).xh(:,end); % update
  tracks(i).P(:,:,end+1) = F*tracks(i).P(:,:,end)*F' + Q; 
  [yEst,Hpre] = Hjac(tracks(i).xh(:,end),exyz);
  tmp = cat(3,Hpre(1:nE,:),Hpre(nE+1:end,:));
  H = cat(2,permute(tmp,[3,2,1]),zeros(2,3,nE)); 
  yExp(:,:,:,i) = permute([yEst(1:nE),yEst(nE+1:end)],[2,3,1]); 
  pExp(:,:,:,i) = mm3d(mm3d(H,tracks(i).P(:,:,end)),permute(H,[2,1,3,4]));
end

% get m distance for detections to tracks
% this can't be parallel: different number of detections per emitter
for i = 1:nE 
  S = det(i).R + pExp(:,:,i,:); % all residual covariances
  res = permute(det(i).y,[2,3,1]) - yExp(:,:,i,:);
  M2 = mm3d(permute(res,[2,1,3,4]),ml3d(S,res));
  det(i).M2 = permute(M2,[3,4,1,2]);
  det(i).res = res;
  det(i).S = S;
end


%% association laws:
% hard association for now (joint probabilistic is warranted)
% track can not take multiple detections from an emitter 
% detection can not go to multiple tracks
pile(nTr).r = [];
pile(nTr).d = [];
pile(nTr).S = [];

for i = 1:nE 
  M2 = det(i).M2;
  pass = M2 < thrBistat & M2 == min(M2) & M2 == min(M2,[],2);
  
  dete = 1:size(pass,1);
  tra = 1:size(pass,2);
  
  [TRA,DETE] = meshgrid(tra,dete);
  tr = TRA(pass);
  de = DETE(pass);
  
  for j = 1:length(tr)
    pile(tr(j)).r(end+1) = det(i).y(de(j),1);
    pile(tr(j)).d(end+1) = det(i).y(de(j),2);
    pile(tr(j)).S(end+1,:) = det(i).s;
  end
  det(i).used(de) = true;
end

% disp('chicken')


%% now get best solution for new measurements 
%update track too
H = [eye(6),zeros(6,3)];
for i = 1:nTr
  
  % add closed form solutions for only 3 associations
  
  
  if length(pile(i).r) < 3
    tracks(i).misses = tracks(i).misses + 1;
    continue
  end
  tracks(i).misses = 0;
  
  [x,P] = itSB(  tracks(i).xh(1:6,end),...
    pile(i).r',pile(i).d',pile(i).S,rErr,dErr);
  P = P*1.5;
  
  res = x - H*tracks(i).xh(:,end);
  S = P + H*tracks(i).P(:,:,end)*H';
  
  K = tracks(i).P(:,:,end)*H'*S^-1;
  tracks(i).xh(:,end) = tracks(i).xh(:,end) + K*res;
  tracks(i).P(:,:,end) = (eye(9)-K*H)*tracks(i).P(:,:,end);
  
end

% remove used detections
for i = 1:nE  
  det(i).R = det(i).R(:,:,~det(i).used);
  det(i).y = det(i).y(~det(i).used,:);
end


% remove bad tracks
keep = true(1,nTr);
for i = 1:nTr
  x = abs(tracks(i).xh(:,end));
  sig = sqrt(abs(diag(tracks(i).P(:,:,end))));
  
  m1 = any(x(4:5)>1000);
  m2 = any(sig(4:5) > 300);
  m3 = tracks(i).misses > 2;
  if m3 || m1 || m2
    keep(i) = false;
  end
  
end
tracks = tracks(keep);

%% end of used code

  % create struct of m2
%   
%   for k = 1:nE % consider every emitter
%     H = [Hpre([k,k+nE],:),zeros(2,3)];
%     yp = yEst([k,k+nE]); % y projected 
%     S = det(k).R + H*M*H'; % residual covariance
%     res = permute(det(k).y,[2,3,1]) - yp;
%     M2 = mm3d(mm3d(permute(res,[2,1,3]),inv2d(S)),res);
%     
%     take = squeeze(M2) < thrBistat;
%     
%     if sum(take) > 1  
%       [~,take] = min(squeeze(M2));
%     end
%     
%     if ~any(take)
%       continue
%     end
%     
%     det(k).used(take) = true;
%     
%     K = M*H'*S(:,:,take)^-1;
%     xp = xp + K*res(:,:,take);
%     M = (eye(9) - K*H)*M;
%   end
% end
%   tracks(j).xh(:,end+1) = xp;
%   tracks(j).P(:,:,end+1) = M;
%   
%   
%   % condition to kill a track
%   
%   tmp = sqrt(abs(diag(M)));
%   if sum(tmp(4:5)) > 1000
%     if isempty(TRACKS)
%       TRACKS = tracks(j);
%     else
%       TRACKS(end+1) = tracks(j);
%     end
%     remove(j) = true;
%   end
%   
% end % loop over tracks
% tracks(remove) = [];
% 
% % return remaining detections for track initializer







