

function [x,P] = itSB(x0,r,d,S,rangeErr,dopErr)

% iteratively minimizes bistatic equations, returns solutions that pass
% threshold test 
% 
% XYZ = permute(xyz,[2,1,3]);
% XYZvel = getVel(XYZ,S,d);
% x0 = permute([XYZ,XYZvel],[2,1,3]);
yData = [r/rangeErr;d/dopErr];

x = x0;

for i = 1:5 % gauss newton iterations 
  [yEst,jac] = Hjac(x,S,rangeErr,dopErr);
  res = yData - yEst;
  jacp = permute(jac,[2,1,3]);
  
%   x = x +mm3d(jacp,res)*1e-3; % this would be gradient descent? 
  
  x = x+ mm3d(  ml3d(mm3d(jacp,jac)    ,jacp   ),res);
end


[~,jac] = Hjac(x,S); % last evaluation of jacobian in non-normalized
jacp = permute(jac,[2,1,3]);


% p = .6;
nE = size(S,1);
% thresh = chi2inv(p,nE*2);
R = diag([ones(1,nE)*rangeErr,ones(1,nE)*dopErr].^2);
Si = diag(1./diag(R));

% res = [r;d] - yEst;
% m2 = mm3d(permute(res,[2,1,3]),mm3d(Si,res));

% idx = m2 < thresh;
tmp  =   gpuArray(mm3d(mm3d(jacp,Si),jac));
P = gather(pagefun(@inv,tmp));
% P = CRB(:,:,idx);

%   y = [r./varargin{1};d./varargin{2}];
%   jac = [elipNorm./varargin{1},zeros(size(elipNorm)); ...
%     dDdX./varargin{2},elipNorm./varargin{2}];

%% usefull equations for posterity 
% Pres = R - jac0*CRB*jac0';






