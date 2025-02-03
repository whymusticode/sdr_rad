

function [x,P,idx,m2] = iterativeSolveBistat(xyz,r,d,S,rangeErr,dopErr)

% iteratively minimizes bistatic equations, returns solutions that pass
% threshold test

% xyz(3,:,xyz(3,:,:)<0) = -xyz(3,:,xyz(3,:,:)<0);

XYZ = permute(xyz,[2,1,3]);
XYZvel = getVel(XYZ,S,d);
x0 = permute([XYZ,XYZvel],[2,1,3]);
yData = [r/rangeErr;d/dopErr];

x = x0;

idx = x(3,:,:) < 0;
x(3,:,idx) = -x(3,:,idx);

i = 0;
maxSteps = 5;
notConverged = true(size(x,3),1,'gpuArray');
while i < maxSteps && any(notConverged) % gauss newton iterations
  i = i+1;
  
%   x = gnStepLoc(x,S,rangeErr,dopErr);
  [yEst,jac] = Hjac(x(:,:,notConverged),S(:,:,notConverged),rangeErr,dopErr);
  res = yData(:,:,notConverged) - yEst;
  jacp = permute(jac,[2,1,3]);
  step = mm3d(  pagefun(@mldivide,mm3d(jacp,jac)    ,jacp   ),res);
  
  
  x(:,:,notConverged) = x(:,:,notConverged) + step;
  crt = squeeze(sum(abs(step./[ones(3,1)*rangeErr;ones(3,1)*dopErr]),1));
  notConverged(notConverged==1) = crt>.1;
  %   x(3,:,17)
end
  %   step = mm3d(jacp,res)*1e-3; % gradient descent


[yEst,jac] = Hjac(x,S); % last evaluation of jacobian in non-normalized
jacp = permute(jac,[2,1,3]);


p = .7;
nE = size(S,1);
thresh = chi2inv(p,nE*2);
R = diag([ones(1,nE)*rangeErr,ones(1,nE)*dopErr].^2);
Si = diag(1./diag(R));

res = [r;d] - yEst;
m2 = mm3d(permute(res,[2,1,3]),mm3d(Si,res));

idx = m2 < thresh;

CRB = gather(pagefun(@inv,gpuArray(mm3d(mm3d(jacp,Si),jac))));
P = CRB(:,:,idx);
x = x(:,:,idx);

%   y = [r./varargin{1};d./varargin{2}];
%   jac = [elipNorm./varargin{1},zeros(size(elipNorm)); ...
%     dDdX./varargin{2},elipNorm./varargin{2}];

%% usefull equations for posterity
% Pres = R - jac0*CRB*jac0';






