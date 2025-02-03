clear

nE = 5; % number of emitters
nT = 1; % number of targets in scene
rangeErr = 100; % [m] rms range error
dopErr = 5; % [m/s] rms doppler error

n = 1e3; % monte carlos

%% generate truth and measurements
exyz = [(rand(nE,2))*1e5,300+100*randn(nE,1)]; % emitter locations
txyz = [(rand(nT,2))*1e5,5e3*ones(nT,1)+10e3*rand(nT,1)]; % targets
vxyz = randn(nT,3)*200; % target velocities

posVel = permute([txyz,vxyz],[2,3,1]);
S = repmat(exyz,[1,1,nT]);
[yTrue,jac0] = Hjac(posVel,S); % truth and jacobian at truth 

y = yTrue + cat(1,randn(nE,1,n)*rangeErr, randn(nE,1,n)*dopErr);
r = y(1:nE,:,:);
d = y(nE+1:end,:,:);

%% initialize with approximation
tic
est = TDOA2XYZgpu(gpuArray([repmat(S,[1,1,n]),r]));
XYZ = permute(est(1:3,:,:),[2,1,3]);
XYZvel = getVel(XYZ,S,d);
pV = permute([XYZ,XYZvel],[2,1,3]); % positionVelocity initialization

%% Gauss Newton
yData = [r/rangeErr;d/dopErr]; % normalize
% for i = 1:5 % grad descent 
%   [yEst,jac] = Hjac(pV,S,rangeErr,dopErr);
%   res = yData - yEst;
%   jacp = permute(jac,[2,1,3]);
%   
%   pV = pV +mm3d(jacp,res)*1e-2;
%   
%   pV = pV + mm3d(  pagefun(@mldivide,mm3d(jacp,jac)    ,jacp   ),res)*i/5;
% end

for i = 1:3 % gauss newton iterations 
  [yEst,jac] = Hjac(pV,S,rangeErr,dopErr);
  res = yData - yEst;
  jacp = permute(jac,[2,1,3]);
  
%   pV = pV +mm3d(jacp,res)*1e-3;
  
  pV = pV + mm3d(  pagefun(@mldivide,mm3d(jacp,jac)    ,jacp   ),res);
end

toc



%% determine acceptance criteria 
p = .6;
thresh = chi2inv(p,nE*2);
R = diag([ones(1,nE)*rangeErr,ones(1,nE)*dopErr].^2);
CRB = (jac0'*R^-1*jac0)^-1;
Pres = R - jac0*CRB*jac0';
[yEst,~] = Hjac(pV,S);
res = [r;d] - yEst;

Si = diag(1./diag(R));
% Smonte1 = cov(permute(res,[3,1,2]))^-1;
m2 = mm3d(permute(res,[2,1,3]),mm3d(Si,res));


% tmp = pagefun(@mldivide,gpuArray(Pres),res);
% m2 = pagefun(@mtimes,permute(res,[2,1,3]),tmp); 

disp(sum(m2 < thresh)/n)


% for ii = 1:gpuDeviceCount
%     g = gpuDevice(ii);
%     fprintf(1,'Device %i has ComputeCapability %s \n', ...
%             g.Index,g.ComputeCapability)
% end

% [V,D] = eig(Pres2);
% sqrt(sort(diag(D)))

%% plot results
xyzDet = permute(pV,[3,1,2]);
idxNew = m2 < thresh;

figure(20)
hold off
errorEllipse(CRB(1:2,1:2)/1e6,txyz(:,1:2)/1e3)
hold on
plot(squeeze(xyzDet(idxNew,1))/1e3,squeeze(xyzDet(idxNew,2))/1e3,'go')
plot(squeeze(xyzDet(~idxNew,1))/1e3,squeeze(xyzDet(~idxNew,2))/1e3,'ro')

font
legend('Cramer Rao Bound','Monte Carlos','Location','EastOutside')
xlabel('East (km)')
ylabel('North (km)')
axis equal
% figure(20)
% hold off
% % plot(exyz(:,1)/1e3,exyz(:,2)/1e3,'d')
% % % plot(txyz(:,1)/1e3,txyz(:,2)/1e3,'+')
% hold on
% plot(squeeze(xyzDet(:,1))/1e3,squeeze(xyzDet(:,2))/1e3,'o')
% font
% legend('Targets','Estimates','Location','EastOutside')
% xlabel('East (km)')
% ylabel('North (km)')



% plot(exyz(:,1)/1e3,exyz(:,2)/1e3,'d')
% plot(txyz(:,1)/1e3,txyz(:,2)/1e3,'+')

% figure(3)
% histogram(squeeze(m2))


%% for TDOA2XYZgpu, if estimate of norm(x) is far from the estimate of x,
%% the overall estimate is likely trash.
% figure(19)
% plot(1:nT,sqrt(sum(xyzDet(:,1:3).^2,2)),1:nT,xyzDet(:,4))

% step 1 TSE
% h = -sum([S.^2,-r.^2],1)/2;
% S2 = [S,-r];
%
% Sp = permute(S2,[2,1,3]);
% Q = repmat(diag(ones(nE,1)*rangeErr^2),[1,1,nT]);
% thetaHat = Sp*

% step 2 TSE
% P = Ppre * rangeErr^2;
% G = [1,0,0;0,1,0;0,0,1;1,1,1];
% D = zeros(4,4,nT,'gpuArray');
% for i = 1:4
%   D(i,i,:) = est(i,:,:);
% end
% O = 4*mm3d(mm3d(D,P),D);
%
% Gp = repmat(G',[1,1,nT]);
% G = repmat(G,[1,1,nT]);
% Oi = pagefun(@inv,O);
% PP = pagefun(@inv,mm3d(mm3d(Gp,Oi),G));
% TH = est.^2;
% x = sqrt(mm3d(mm3d(mm3d(PP,Gp),Oi),TH));
% xyzDet = permute(x,[3,1,2]);


%   hes = mm3d(jacp,jac); % hessian
%   hesDiag = hes.*repmat(eye(size(hes,1)),[1,1,size(hes,3)]);
%   pV = pV + mm3d(mm3d(pagefun(@inv,hesDiag*lam + hes),jacp),res);






% (eye(8)*1e-10 + Pres)^-1
% bb = - jac*CRB*jac';
% aa = R;
% aai = diag(1./diag(R));
% norm(aai*bb)
% g = trace(bb*aai);
% Si = aai - aai*bb*aai + aai*bb*aai*bb*aai;
% Si*Pres
% Si*

% Si = aai - aai*bb*(eye(8) + aai*bb)^-1*aai;

% return
% [V,D] = eig(Pres);




% S = ones(3)+eye(3);
% [V,D] = eig(S);
% tmp = sqrt(diag(D));
% V2 = V.*tmp';
% D2i = diag(1./diag(D).^2);
% Si = V2'^-1*D2i*V2^-1;

% return

% [V,D] = cdf2rdf(V,D);


% Vi = real(V^-1);
% V = real(V);
% D = real(D);
% Di = diag(1./diag(D));
% Si = V*Di*V^-1;


% res(:,:,i)'*(Pres\res(:,:,i))
%
% S = ones(3)+eye(3);
% [V,D] = eig(S);
% tmp = sqrt(diag(D));
% V*D*V'
% V2 = V.*tmp';
% V2*V2'
% V2'^-1*V2^-1



% V2'^-1*V2^-1


% x = [1;2;3];
% b = [-.5;.5;1.5];
% S^-1*x
% S\x
% (Pres\res(:,:,i))



% Si = Pres^-1;
% Si = R^-1;