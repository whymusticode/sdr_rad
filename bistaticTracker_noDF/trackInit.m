function [tracks,trNew] = trackInit(det,trPar,tracks)
% det(1).y is all the detections for 1st emitter
% (range and doppler measured)

% if too many detections get passed, this function will be a computational
% disaster, as it scales with the #of combinations of detections

trNew = struct([]);

rangeErr = trPar.rErr;
dopErr = trPar.dErr;
nReq = 4; % number of emitter detections required [4 or 5]

[S,r,d,bigidx,nPerEmit] = makeCombos(det,nReq); %#ok<ASGLU>

if isempty(S);return;end

xyz = coldStart([S,r]); % initialization solves for 4 variables
xyz = xyz(1:3,:,:); % 4th variable is useless
% might be wise to neck down here in some cases


%% Iterative solver
[x,P,idx,m2] = iterativeSolveBistat(xyz,r,d,S,rangeErr,dopErr);
x = gather(x);
P = gather(P);
idx = squeeze(gather(idx));
m2 = gather(m2);
% I think ever different integer in idx is a specific detection

%% only keep best solution for each detection
idxMat = bigidx(idx,:);
% x = x(:,:,idx);
% P = P(:,:,idx);
idxAll = [];

m3 = m2(idx);
while size(idxMat,1) > 0
  [~,loc] = min(m3);
  tmp = false;
  for j = 1:nReq
    tmp = any(idxMat == idxMat(loc,j),2) | tmp;
  end
  
  idxAll = [idxAll,loc]; %#ok<AGROW>
  m3 = m3(~tmp);
  idxMat = idxMat(~tmp,:);
end


x = x(:,:,idxAll);
P = P(:,:,idxAll);
% disp([numel(idx),sum(idx),numel(idxAll)])

if isempty(x);return;end

for i = 1:size(x,3)
  
% %   throw out anything funny at initialization
  xxx = abs(x(:,:,i));
  sig = sqrt(abs(diag(P(:,:,i))));
  m1 = any(xxx(4:5)>700);
  m2 = any(sig(4:5) > 100);
  if m1 || m2
    continue
  end
  
  
  tracks(end+1).xh = [x(:,:,i);zeros(3,1)]; %#ok<AGROW>
  tracks(end).P = P(:,:,i);
  tracks(end).P(7:9,7:9) = eye(3)*trPar.maxAcc^2;
  tracks(end).misses = 0;
  tracks(end).color = rand(3,1);
  
  if isempty(trNew)
    trNew = tracks(end);
  else
  trNew(end+1) = tracks(end);
  end
end




end





function [S,r,d,bigidx,nPerEmit] = makeCombos(det,nReq)
% create all combinations of nReq (4) detections

% I bet there is a way to do this more efficiently, is pretty slow now

maxN = 2e5; % max number of combos to consider (4gb gpu ram for 2e6)

nE = length(det);
nPerEmit = zeros(1,nE);
for i = 1:nE % get # of detections for each emitter
  nPerEmit(i) = size(det(i).y,1);
end

% count number of combinations
comboIdx = nchoosek(1:nE,nReq); % choose nReq emitters
nCombs = sum(prod(nPerEmit(comboIdx),2)); % total number of combinations
if nCombs > maxN
  disp(nCombs)
  error('too many combinations')
end

if nCombs < 1; S = []; r = []; d =[]; bigidx = [];return;end


cumPer = cumsum(nPerEmit)-nPerEmit;
%% create combinations (This is slow, and should be improved)
% tic
comboIdx = nchoosek(1:nE,nReq); % choose nReq emitters
r = []; % measured bistatic range matrix
S = []; % emitter position matrix
d = []; % doppler matrix
bigidx = []; % stores emitter number and detection id for that emitter
for i = 1:size(comboIdx,1)
  IDs = comboIdx(i,:);
  sizes = nPerEmit(IDs);
  switch nReq
    case 5
      [X,Y,Z,W,WW] = ndgrid(1:sizes(1),1:sizes(2),1:sizes(3),...
        1:sizes(4),1:sizes(5));
      XYZW = [X(:),Y(:),Z(:),W(:),WW(:)];
    case 4
      [X,Y,Z,W] = ndgrid(1:sizes(1),1:sizes(2),1:sizes(3),...
        1:sizes(4));
      XYZW = [X(:),Y(:),Z(:),W(:)];
    otherwise
      error('nReq can only be 4 or 5')
  end
  n = size(XYZW,1);
  rtmp = zeros(4,1,n,'single','gpuArray');%,'single','gpuArray'
  dtmp = rtmp;
  sPre = zeros(4,3,1,'single','gpuArray');%,'single','gpuArray'
  for j = 1:4
    rtmp(j,:,:) = permute(det(IDs(j)).y(XYZW(:,j),1),[2,3,1]);
    dtmp(j,:,:) = permute(det(IDs(j)).y(XYZW(:,j),2),[2,3,1]);
    sPre(j,:) = det(IDs(j)).s;
  end
  stmp = repmat(sPre,[1,1,n]);
  r = cat(3,r,rtmp); % r is 4 x 1 x nCombs
  d = cat(3,d,dtmp);
  S = cat(3,S,stmp); % S is 4 x 3 x nCombs
  idxtmp = cumPer(IDs)+XYZW;
  bigidx = cat(1,bigidx,idxtmp);
end
end





%% end used code


% for i = 1:4 % gauss - newton
%   [yEst,jac] = Hjac(pV,S,rangeErr,dopErr);
%   res = yData-yEst;
%   jacp = permute(jac,[2,1,3]);
%   pV = pV + mm3d(  pagefun(@mldivide,mm3d(jacp,jac), jacp ),res);
% end
%
% errF = squeeze(sum(res.^2,1))/(2*nReq);
% finDx = errF < 2;
% xyzDet = permute(pV(:,:,finDx),[3,2,1]);
% bigidx = bigidx(finDx,:);
% errF = errF(finDx);
%
%
%
% if isempty(pile);return;end
%
% % format the initializations
% for i = 1:length(pile)
%   pile(i).xh = [pile(i).xyz';zeros(3,1)];
%   pile(i).P(7:9,7:9) = eye(3)*trPar.maxAcc^2;
% end
%
%
%
% if isempty(tracks)
%   tracks = pile;
% else
%   tracks(end+1:end+length(pile)) = pile;
% end




%% neck down
% max2Iterate = 1e4; % most number of combinations allowed through next step
% if nCombs > max2Iterate
%
%   pos = permute(xyz,[2,1,3]);
%   xms = pos - S;
%   smxn = sqrt(sum(xms.^2,2));
%   xn = sqrt(sum(pos.^2,2));
%   err = smxn + xn - r;
%   err2 = squeeze(sum(err.^2,1));
%   errS = sort(err2); % obv really slow
%   idx = err2 <= errS(max2Iterate);
%
%   S = S(:,:,idx);
%   xyz = xyz(:,:,idx);
%   r = r(:,:,idx);
%   d = d(:,:,idx);
%   bigidx = bigidx(idx,:);
% end


%% EXTENDED ASSOCIATION
% this was a really complicated way of gaurunteeing no measurements used
% twice, and combining
%
%
% cumPerC = cumsum(nPerEmit);
% cumPer = cumPerC - nPerEmit;
%
% bigDetMat = [];
% for i = 1:length(det)
%   bigDetMat = [bigDetMat;det(i).y]; %#ok<AGROW>
% end
% nE = length(det);
% Sall = zeros(nE,3);
% for i = 1:nE % put into realistic format with noise
%   Sall(i,:) = det(i).s;
% end
%
% pile = struct([]);
% testDx = 1;
% opti = optimoptions('lsqcurvefit',...
%   'SpecifyObjectiveGradient',true); %#ok<NASGU>
% func = @(x,xd)Hjac(x,xd,rangeErr,dopErr);
% while ~ isempty(bigidx) && testDx < size(bigidx,1)
%   myset = any(permute(bigidx(testDx,:),[1,3,2]) == bigidx,3);
%   assocIdx = sum(myset,2) == 3;
%
%   if sum(assocIdx) == 0
%     testDx = testDx + 1;
%     continue
%   end
%
%   potAssoc = unique(bigidx(~myset & assocIdx))'; % potential associations
%   alle = sort([bigidx(testDx,:),potAssoc]);
%   matchDetEm = alle > cumPer' & alle <= cumPerC';%[emitter x detect]
%   assocPerEm = sum(matchDetEm,2);
%
%   full = {};
%   for i = 1:nE
%     tmp = alle(matchDetEm(i,:));
%     if length(tmp) > 1
%       full{end+1} = tmp; %#ok<AGROW>
%     end
%   end
%   if isempty(full)
%     potentAddOn = [];
%   else
%     n = numel(full); % number of vectors
%     combs = cell(1,n); % pre-define to generate comma-separated list
%     [combs{end:-1:1}] = ndgrid(full{end:-1:1}); % the reverse order in these two
%     % comma-separated lists is needed to produce the rows of the result matrix
%     combs = cat(n+1, combs{:}); %concat the n n-dim arrays along dimension n+1
%     potentAddOn = reshape(combs,[],n); %reshape to obtain desired matrix
%   end
%
%   alle = alle(~any(matchDetEm&(assocPerEm>1),1));
%
%   xdata = [];
%   ydata = [];
%   init = double(gather(permute(xyzDet(testDx,:,:),[3,1,2]))); %#ok<NASGU>
%   %   init = permute(xyzDet(testDx,:,:),[3,1,2]);
%   if isempty(potentAddOn)
%     %     init = double(gather(permute(xyzDet(testDx,:,:),[3,1,2]))); %#ok<NASGU>
%     allen = alle;
%     ydata = [bigDetMat(allen,1)./rangeErr;bigDetMat(allen,2)./dopErr];
%     xdata = Sall(assocPerEm>=1,:);
%     T = evalc(['[sol,resnorm] = lsqcurvefit(func,init,xdata',...
%       ',ydata,[],[],opti);']); %#ok<NASGU>
%   else % look over potential add ons to a track
%
%     resnorm2 = [];
%     sol2 = {};
%     for i = 1:size(potentAddOn,1)
%       allen = [alle,potentAddOn(i,:)];
%       matchDetEm = allen > cumPer' & allen <= cumPerC';%[emitter x detect]
%       assocPerEm2 = sum(matchDetEm,2);
%       ydata = [bigDetMat(allen,1)./rangeErr;bigDetMat(allen,2)./dopErr];
%       xdata(:,:,i) = Sall(assocPerEm2>=1,:); %#ok<AGROW>
%       T = evalc(['[sol2{i},resnorm2(i)] = lsqcurvefit(',...
%         'func,init,xdata(:,:,i),ydata,[],[],opti);']); %#ok<NASGU>
%     end
%     [resnorm,loc] = min(resnorm2);
%     sol = sol2{loc};
%     alle = [alle,potentAddOn(loc,:)]; %#ok<AGROW>
%     xdata = xdata(:,:,loc);
%   end
%   if resnorm/length(ydata) > 3 % solution can not be conflated with others
%     testDx = testDx + 1; % try next solution
%     continue
%   end
%   % we now have a new high confidence solution
%   pile(end+1).det = alle; %#ok<AGROW> % indicies
%   pile(end).xyz = sol';
%   CRB = getCrb(func,sol,xdata,ones(size(ydata)));
%   pile(end).P = CRB; % iterative solver should perform at CRB
%
%   % remove anything with a detection used by the new solution
%   idx2 = ~any(any(permute(alle,[1,3,2]) == bigidx,3),2);
%   testDx = testDx - sum(~idx2(1:testDx)) + 1;
%   bigidx = bigidx(idx2,:);
%   xyzDet = xyzDet(idx2,:,:);
%   errF = errF(idx2,:,:);
% end
%
% % add any combination of only 4 that might be left
% for i = 1:size(bigidx,1)
%   alle = bigidx(i,:);
%   matchDetEm = alle > cumPer' & alle <= cumPerC';%[emitter x detect]
%   assocPerEm = sum(matchDetEm,2);
%   pile(end+1).det = bigidx(i,:); %#ok<AGROW> % indicies
%   pile(end).xyz = double(gather(permute(xyzDet(i,:,:),[3,2,1])))';
%   xdata = Sall(assocPerEm>=1,:);
%   ydata = [bigDetMat(bigidx(i,:),1)./rangeErr;...
%     bigDetMat(bigidx(i,:),2)./dopErr];
%   CRB = getCrb(func,pile(end).xyz',xdata,ones(size(ydata)));
%   pile(end).P = CRB; % iterative solver should perform at CRB
% end















%     for i = 1:size(potentAddOn,1)
%       allen = [alle,potentAddOn(i,:)];
%       matchDetEm = allen > cumPer' & allen <= cumPerC';%[emitter x detect]
%       assocPerEm2 = sum(matchDetEm,2);
%       ydata(:,:,i) = [bigDetMat(allen,1)./rangeErr;...
%         bigDetMat(allen,2)./dopErr];
%       xdata(:,:,i) = Sall(assocPerEm2>=1,:); %#ok<AGROW>
%     end
%     lam = .001; % damping rate
%     pV = init;
%     for i = 1:5
%       [yEst,jac] = Hjac(pV,xdata,rangeErr,dopErr);
%       res = ydata-yEst;
%
%       % Levenburg-Marquardt
%       jacp = permute(jac,[2,1,3]);
%       hes = mm3d(jacp,jac); % hessian
%       hesDiag = hes.*repmat(eye(size(hes,1)),[1,1,size(hes,3)]);
%       upd = mm3d(mm3d(pagefun(@inv,hesDiag*lam + hes),jacp),res);
%       pV = pV + upd;
%     end
%
%     resnorm2 = squeeze(sum(res.^2,1));
%     [resnorm,loc] = min(resnorm2);
%     sol = pV(:,:,loc);


%% evaluate initialization for all combinations
% this section is based on "Two Methods for Target Localization in
% Multistatic Passive Radar" by MATEUSZ MALANOWSKI and KRZYSZTOF KULPA
% z = sum([S.*S,-r.*r],2)/2;
% Sp = permute(S,[2,1,3]);
% S2 = mm3d(Sp,S);
% s2i = pagefun(@inv,S2);
% Sst = mm3d(s2i,Sp); % s star
% a = mm3d(Sst,z);
% b = mm3d(Sst,r);
%
% if 0 %spherical intersection method of getting RtHat:
%   m = size(S,1); %#ok<UNRCH> % number of transmitters
%   rp = permute(r,[2,1,3]);
%   T = eye(m,'gpuArray') - mm3d(S,Sst);
%   rtt = mm3d(rp,T);
%   RtHat = - mm3d(rtt,z)./mm3d(rtt,r); % SX method
%   idx = true(size(S,3),1);
% else %"spherical interpolation" method
%   ap = permute(a,[2,1,3]);
%   bp = permute(b,[2,1,3]);
%   btb1 = (mm3d(bp,b)-1);
%   determinant = 4*mm3d(ap,b).^2 - 4*btb1.*mm3d(ap,a);
%   determinant(determinant<0) = 0;
%   part1 = -2*mm3d(ap,b);
%   sqrtDet = sqrt(determinant);
%   RtHat = (part1 + sqrtDet)./(2*btb1);
%   RtHat2 = (part1 - sqrtDet)./(2*btb1);
%   RtHat(RtHat<0) = RtHat2(RtHat<0); % ??
% end
% xyz = a + b.*RtHat;


%% other shity code
% return
% y = Hbistat(XYZ,XYZvel,S); % gets what we should have measured

% figure(20)
% hold on
% plot(squeeze(XYZ(:,1,:)),squeeze(XYZ(:,2,:)),'o')

% posErr = (r - y(:,1,:))/rangeErr;
% velErr = (d- y(:,2,:))/dopErr;


% disp(permute(posErr,[1,3,2]))
% disp(permute(velErr,[1,3,2]))

% err = sum((permute(posErr,[1,3,2])).^2,1) + ...
%   sum((permute(velErr,[1,3,2])).^2,1);

% [~,wow] = sort(err)


% initPos = XYZ(:,:,subIdx);
% initVel = XYZvel(:,:,subIdx);
% func = @(X,S)Hbistat(X(1,:,:),X(2,:,:),S,rangeErr,dopErr);
% xData = double(gather(S(:,:,subIdx)));
% bigInit = double(gather(cat(1,initPos,initVel)));
% yData = double(gather(cat(2,r(:,:,subIdx)/rangeErr,d(:,:,subIdx)/dopErr)));
% func(bigInit,xData); % for testing purposes
%
% errPP = err;
% thresh = 1; % number of standard deviations for "solution"
% remDet = [];
%
% tic
% while 1
%   if toc > 5
%     break
%   end
%
%   [~,loc] = min(errPP);
%
%   % iterative solver, use
%   init = double(gather([XYZ(:,:,loc);XYZvel(:,:,loc)]));
%   xData = double(gather(S(:,:,loc)));
%   yData = double(gather([r(:,:,loc)/rangeErr,d(:,:,loc)/dopErr]));
%   T = evalc('[sol,err2] = lsqcurvefit(func,init,xData,yData);');



%   if err2 < nReq*2*thresh % angels rejoice, we have a solution
%     remDet(:,:,end+1) = bigidx(:,:,loc);
%     idxF = removeDetects(bigidx,remDet(:,:,end));
%     errPP(idxF) = 1e20;
%   else
%     errPP(loc) = 1e20; % this combo will be removed from consideration
%
%   end
%
%
%
% end




%   eval(:,i) = totErr([1;16]);

%   [valdc,idc] = sort(totErr);
%   idc(1:2)

%   histogram(squeeze(sum(res.^2,1)))
%   pause(1)

%% other fast iterative solvers
% gauss - newton (sucks)
% jacp = permute(jac,[2,1,3]);
% pV = pV - mm3d(mm3d(pagefun(@inv,mm3d(jacp,jac)),jacp),res);

% grad descent (slow)
%   pV = pV + mm3d(permute(jac,[2,1,3]),res)*1e1;% grad descent


%% doing it this way supports any value for nReq, but it's slow
%   in = ''; out = '';varName = 'ABCDEFGHI';cond = '';
%   for j = 1:nReq
%     in = [in,',1:sizes(',num2str(j),')'];
%     out = [out,',',varName(j)];
%     cond = [cond,',',varName(j),'(:)'];
%   end
%   T = evalc(['[',out(2:end),'] = ndgrid(',in(2:end),')']);
%   T = evalc(['XYZW = [',cond(2:end),']']);

%% this is just meant to verify the derivation of jacobian
% [yEst1,jac] = Hjac(pV,S,rangeErr,dopErr);
% pv = pV;
% dx = 1;
% nIdx = 10;
% xyidx = 6;
% pv(xyidx,1,nIdx) = pv(xyidx,1,nIdx) + dx;
% [yEst2,~] = Hjac(pv,S,rangeErr,dopErr);
% disp([(yEst2(:,:,nIdx)-yEst1(:,:,nIdx))/dx,jac(:,xyidx,nIdx)])

%% can we leverage matlab's algorithm for better performance?
% options = optimoptions('lsqcurvefit','SpecifyObjectiveGradient',true);
% myFun = @(x,xd)Hjac(x,xd,rangeErr,dopErr);
% xData = double(gather(S));
% yData = double(gather(y));
% init = double(gather(pV));
% sol = lsqcurvefit(myFun,init,xData,yData,[],[],options);


% feat = sum([abs(permute(posErr,[1,3,2]));...
%   abs(permute(velErr,[1,3,2]))],1);


% figure(46)
% plot(pErr,vErr,'.')
% grid on

% histogram(log(sum(vErr,1)))


% xyzDet = [permute(XYZ(:,:,subIdx),[3,2,1]),...
%   permute(XYZvel(:,:,subIdx),[3,2,1])];


% subIdx1 = permute(sum(posErr.^2)/4,[3,1,2]) < inf;
% subIdx2 = permute(sum(velErr.^2)/4,[3,1,2]) < inf;
% subIdx = subIdx1|subIdx2;
% subIdx = permute(sum(posErr.^2 + velErr.^2),[3,1,2]) < NNN*8;
% sum(velErr.^2)/4



% sum(subIdx)
% % tmp = velErr(:,:,subIdx);
% % min(tmp(:))
%
% y0 = y(:,:,subIdx);
%
% hold off
% plot(txyz(:,1),txyz(:,2),'+')
% hold on
% plot(squeeze(XYZ(1,1,subIdx)),squeeze(XYZ(1,2,subIdx)),'o')
%
% hold off
% plot(vxyz(:,1),vxyz(:,2),'+')
% hold on
% plot(squeeze(XYZvel(1,1,subIdx)),squeeze(XYZvel(1,2,subIdx)),'o')





% subidxtmppre = bigidx(:,:,idx);
% subidxtmp = subidxtmppre(:,:,subIdx);
%
% % sum(subIdx)
% cartesianDetections = permute(xyz(:,:,subIdx),[3,1,2]);





%% lastly, perform iterative solution for most accurate results
% nIterate = sum(subIdx);
% func = @(X,S)Hbistat(X(1,:,:),X(2,:,:),S,rangeErr,dopErr);
% xData = double(gather(S(:,:,subIdx)));
% bigInit = double(gather(cat(1,initPos,initVel)));
% yData = double(gather(cat(2,r(:,:,subIdx)/rangeErr,d(:,:,subIdx)/dopErr)));
% % func(bigInit,xData); % for testing purposes
%
% for i = 1:nIterate
% T = evalc(['sol(:,:,i) = lsqcurvefit(func,bigInit(:,:,i),',...
% 'xData(:,:,i),yData(:,:,i));']);
% end




% err = yData - func(sol,xData);


% err2 = sum(err.^2,1);

% figure(46)
% loglog(squeeze(err2(:,1,:)),squeeze(err2(:,2,:)),'.')
% grid on
% newDx = squeeze(sum(err2,2)) < 1;
%
%
% xyzDet = [permute(sol(1,:,newDx),[3,2,1]),...
%   permute(sol(2,:,newDx),[3,2,1])];



%% end of useful code



%% working algorithm
% pErr = log10(sum((permute(posErr,[1,3,2])).^2,1));
% vErr = log10(sum((permute(velErr,[1,3,2])).^2,1));
% sssu = 0;
% yInt = -20;
% while sssu < 200
%   yInt = yInt + .1;
%   subIdx = vErr < -pErr + yInt;
%   sssu = sum(subIdx);
% end
% initPos = XYZ(:,:,subIdx);
% initVel = XYZvel(:,:,subIdx);
% nIterate = sum(subIdx);
% func = @(X,S)Hbistat(X(1,:,:),X(2,:,:),S,rangeErr,dopErr);
% xData = double(gather(S(:,:,subIdx)));
% bigInit = double(gather(cat(1,initPos,initVel)));
% yData = double(gather(cat(2,r(:,:,subIdx)/rangeErr,d(:,:,subIdx)/dopErr)));
% for i = 1:nIterate
% T = evalc(['sol(:,:,i) = lsqcurvefit(func,bigInit(:,:,i),',...
% 'xData(:,:,i),yData(:,:,i));']);
% end
% err = yData - func(sol,xData);
% err2 = sum(err.^2,1);
% newDx = squeeze(sum(err2,2)) < 1;
% xyzDet = [permute(sol(1,:,newDx),[3,2,1]),...
%   permute(sol(2,:,newDx),[3,2,1])];







%   [X,Y,Z,W] = ndgrid(1:sizes(1),1:sizes(2),1:sizes(3),1:sizes(4));
%   XYZW = [X(:),Y(:),Z(:),W(:)];

% for spherical intersection method of getting RtHat:
% m = size(S,1); % number of transmitters
% rp = permute(r,[2,1,3]);
% T = repmat(eye(m,'gpuArray'),[1,1,n]) - mm3d(S,Sst);
% rtt = mm3d(rp,T);
% RtHat = - mm3d(rtt,z)./mm3d(rtt,r); % SX method
