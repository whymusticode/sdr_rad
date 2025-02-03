
% This code proves out the ability to use >3 bistatic detection to solve
% ambiguity problems

% computation time scales with nT ^ nE
clear
addpath(genpath('C:\Users\bestm\Desktop\MATLAB\FILTERING'))

nSteps = 10;
nParticles = 1e3;
ts = .5; % seconds integration time, seconds between detections

nE = 7; % number of emitters
nT = 6; % number of targets in scene
rangeErr = 50; % [m] rms range error
dopErr = 10; % [m/s] rms doppler error
prob = 1; % probability of one detection getting through
%% generate random positions
exyz = [(rand(nE,2)-.5)*1e5,300+100*randn(nE,1)]; % emitter locations
txyz = [(rand(nT,2)-.5)*1e5,5e3*ones(nT,1)+10e3*rand(nT,1)]; % target locations
% txyz = [(rand(nT,2)-.5)*1e3+2e4,10e3*ones(nT,1)]; % close targets
vxyz = [randn(nT,2)*200,zeros(nT,1)]; % target velocities

X = zeros(dims,nSteps);
P = zeros(dims,dims,nSteps);
for j = 1:nSteps
  
  y = Hbistat(permute(txyz+vxyz*ts*j,[3,2,1]),permute(vxyz,[3,2,1]),...
    exyz); % get all possible exact measurements
  
  idx = rand(nT,nE) < prob;
  % idx = [1,1,1;0,1,1;0,0,1]==1;
  
  invis = sum(sum(idx,2)<4);
  
  nPerEmit = sum(idx,1);
  for i = 1:nE % put into realistic format with noise
    det(i).s = exyz(i,:);
    det(i).rrd = permute(y(i,:,idx(:,i)),[3,2,1]) + ...
      randn(nPerEmit(i),2).*[rangeErr,dopErr];
  end
  % throw in a false alarm
%   det(3).rrd = [det(3).rrd;...
%     min(det(3).rrd)+rand(1,2).*(max(det(3).rrd)-min(det(3).rrd))];
%   nPerEmit(3) = nPerEmit(3) + 1;
  % profile on
  
  
  % Filtering loop for bootstrap filter
   SX = ungm_f(SX,k) + gauss_rnd(0,1,size(SX,2));
   W  = gauss_pdf(Y(:,k),ungm_h(SX,k),1);
   ind = resampstr(W);
   SX = SX(:,ind);
   M = mean(SX);
   P = var(SX);
   X(:,j)   = M;
   P(:,:,j) = P;    

  
  
  
end



return

pile = bistatDisambig(det,rangeErr,dopErr);
% profile viewer




% comparison with truth:
cumPerC = cumsum(nPerEmit);
cumPer = cumPerC - nPerEmit;
tmp = cumsum(idx)+cumPer;
truth = struct([]);
for i = 1:nT
  truth(end+1).det = tmp(i,idx(i,:));
  truth(end).xyz = [txyz(i,:),vxyz(i,:)];
end

testCRBass = [];
disp('check, any unmatched cases:')
for i = 1:length(pile)
  foundd = 0;
  for j = 1:length(truth)
    
    tr = sort(truth(j).det);
    pi = sort(pile(i).det);
    if numel(pi) == numel(tr)
      if sum(pi == tr) == numel(tr)
        res = truth(j).xyz - pile(i).xyz;
        disp(res)
        disp(sqrt(diag(pile(i).P))')
        testCRBass = [testCRBass,res./sqrt(diag(pile(i).P))'];
        foundd = 1;
      end
    end
  end
  
  if foundd == 0
    disp(num2str(i))
    disp('this pile was not matched in truth')
  end
end

pid = zeros(length(pile),nE);
for i = 1:length(pile)
  pid(i,:) = [pile(i).det,zeros(1,nE-length(pile(i).det))];
end
tid = zeros(length(truth),nE);
for i = 1:length(truth)
  tid(i,:) = [truth(i).det,zeros(1,nE-length(truth(i).det))];
end

%% plot results
figure(20)
hold off
plot(exyz(:,1)/1e3,exyz(:,2)/1e3,'d')
hold on
plot(txyz(:,1)/1e3,txyz(:,2)/1e3,'+')
% plot(squeeze(xyzDet(:,1))/1e3,squeeze(xyzDet(:,2))/1e3,'o')
estXyz = [];
for i = 1:length(pile)
  estXyz(i,:) = pile(i).xyz;
end
plot(squeeze( estXyz(:,1))/1e3,squeeze( estXyz(:,2))/1e3,'o')
font
legend('Transmitters','Targets','Estimates','Location','EastOutside')
xlabel('East (km)')
ylabel('North (km)')


figure(21)
histogram(testCRBass)
title(num2str(std(testCRBass)))




return











%%
truPV = [txyz,vxyz];
loca = [];
for i = 1:nT
  [~,loca(i)] = min(sum((truPV(i,:)-xyzDet).^2,2));
end
truPV - xyzDet(loca,:)

zz = 0:3;
sum(.8.^zz.*.2.^(8-zz).*[nchoosek(8,0),nchoosek(8,1),...
  nchoosek(8,2),nchoosek(8,3)])



%% this is for plotting
% [xx,yy,matches] = drawContours(det);
%
% [EM,TG] = meshgrid(1:nE,1:nT);
% tgidx = TG(idx);
% emidx = EM(idx);
%
% figure(34)
% hold off
% imagesc(xx,yy,matches)
% hold on
% plot(txyz(:,1)/1e3,txyz(:,2)/1e3,'k+','markerSize',10)
% plot(exyz(:,1)/1e3,exyz(:,2)/1e3,'kd','markerSize',10)
%
% % for i = 1:length(tgidx)
% %   line([txyz(tgidx(i),1)/1e3,exyz(emidx(i),1)/1e3],...
% %     [txyz(tgidx(i),2)/1e3,exyz(emidx(i),2)/1e3],'color','g')
% % end
% line([txyz(tgidx,1)/1e3,exyz(emidx,1)/1e3]',...
%   [txyz(tgidx,2)/1e3,exyz(emidx,2)/1e3]','color','g')
%
% tgidx = TG(~idx);
% emidx = EM(~idx);
% line([txyz(tgidx,1)/1e3,exyz(emidx,1)/1e3]',...
%   [txyz(tgidx,2)/1e3,exyz(emidx,2)/1e3]','color','r')
% legend('Targets','Emitters','detection','no detection')
% xlabel('E (km)')
% ylabel('N (km)')
% axis equal
%
% cm = colormap;
% cm(1,:) = 1;
% colormap(cm)
% font
% return

%% for old estimation:

% err = permute(xyz,[3,1,2]) - txyz;
% sqrt(sum(err(:,1:2).^2,2))

% figure(34)
% hold off
% plot(squeeze(xyz(1,1,:)),squeeze(xyz(2,1,:)),'o')
% hold on
% plot(txyz(:,1),txyz(:,2),'+')
% font

% return

%% Use iterative solver to evaluate solutions more accurately
% midx = (xyz(3,1,:) < 0)+1;
%
% sidx = S(:,1:3,idx);
% EXYZ = permute(repmat(sidx,[1,1,2]),[2,1,3]);
% x = cat(3,xyz(1:3,1,:),xyz(1:3,2,:));
%
% dx2  = [dx(:,idx),dx(:,idx)];
% rdoaEst = reshape(f(x,EXYZ),[nE,size(dx2,2)]); % rdoa for
%
% rdoaERR = sum((rdoaEst - dx2).^2);
% [val,sortIdx] = sort(rdoaERR);

%% eliminate up to a certain number of combinations by accuracy so far
% a = 40000;
% try2get = 8*nT; % number of combinations that will be evaluated closer
% good2go = 0;
% count = 0;
% while ~good2go && count < 100
%   count = count + 1;
%   newDx = rdoaERR < rangeError^2*nE*a & squeeze(x(3,:,:)>0)';
%   if abs(sum(newDx) - try2get)/try2get < .5
%     good2go = 1;
%   end
%   if sum(newDx) > try2get
%     a = a/1.5;
%   else
%     a = a*1.5;
%   end
% end
%
%
% if sum(newDx) > 4e2
%   error('too many possible solutions')
% end
%
% tmp = x(:,:,newDx);
% x2 = gather(reshape(tmp,[3,size(tmp,3)]));
% fEXYZ = gather(EXYZ(:,:,newDx));
% finalRD = dx2(:,newDx);
% % matlab least square solver can be shown to perform at CRB, I use that to
% % get most accurate estimates, then make final removal step
% tic
% for i = 1:size(x2,2)
% y = finalRD(:,i)';
%  T = evalc([' [lsqX(:,i),normi(i)] = lsqcurvefit(f,x2(:,i)',...
%    ',fEXYZ(:,:,i),y)']);
% end
% iterTime = toc;
% disp([num2str(sum(newDx)),' combinations evaluated with '...
%   'iterative solver in  ',num2str(iterTime),' seconds']  )
%
%
% detections = lsqX(:,normi<4*rangeError^2*nE)'; % final result




%% extra code

%     xp = x';
%     TT = evalc('[xp,resnorm] = lsqcurvefit(f,xp,[],dx);');
% disp(resnorm)
%     disp('found target.  residual measurement error [e,n,u]:')

%% format of RDOA function for GN minimization
% f = @(tg,~)(sqrt(sum((exyz'-tg).^2))+ sqrt(sum(tg.^2)) - ...
%   sqrt(sum(exyz'.^2)))';

%% code for sharing
% clear
% nE = 5; % number of emitters
% rangeError = 10; %rms range error
%
% exyz = [rand(nE,2)*1e5,300+50*randn(nE,1)]; % emitter locations
% txyz = [rand(1,2)*1e5,rand(1,1)*1e3+2e3]; % target locations
%
% e2s = sqrt(sum(exyz.^2,2)); % range emitter to sensor
% rdoa = @(tg,~) sqrt(sum((exyz'-tg).^2))+ sqrt(sum(tg.^2)) - e2s';
% rdoaMeas = rdoa(txyz') + randn(nE,1)*rangeError;
% xyzEst = TDOA2XYZ([exyz,rdoaMeas+e2s]);

% for i = 1:6
% myz(i,:) = reshape(xyz(1:3,midx(i),i),[1,3]);
% clock(i) = xyz(4,midx(i),i);
% end
% myz - txyz
%
% J = @(tg,exyz,rdoa)sum((sqrt(sum((exyz'-tg).^2))+ sqrt(sum(tg.^2)) - ...
%   sqrt(sum(exyz'.^2))-rdoa).^2);% cost function for grad descent
% JJ = @(tg)J(tg,exyz,rdoa(3,:));
%
% a = 1e-2;
% xx0 = gather(myz(3,:)');
% for i = 1:30
%   gradd = GD(JJ,xx0);
%   xx0 = xx0-gradd*a;
%   disp(JJ(xx0))
% end
%
% JJ(txyz(3,:)')
%
% xx0
% txyz(3,:)-gather(myz')'
%
% -42.9806   43.8514  339.5963

% t3s = sqrt(sum(myz.^2,2)); % range target to sensor
% rdoa2 = zeros(nT,nE,'gpuArray'); % all the range difference measurements
% for i = 1:nE
%   rdoa2(:,i) = sqrt(sum((myz-exyz(i,:)).^2,2)) + t3s - e2s(i);
% end
%
% rdoa2-rdoa
% return
% myz


% rdoaERR2 = reshape(f(reshape(txyz',[3,1,6]),repmat(exyz',[1,1,6])),...
%   [5,6])-dx(:,idx);
