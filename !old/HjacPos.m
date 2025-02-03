function [y,jac] = HjacPos(pos,S,varargin)
% IN:
% posVel  6x1xN   Position and velocity [x;y;z;xd;yd;zd]
% s       Mx3xN   emitter positions
% OUT:
% y - 2Mx1xN   (M emitters * 2 is range then doppler,1, N samples)
% jac - 2Mx6xN (6 is position then velocity)

pos = permute(pos(1:3,:,:),[2,1,3]);

xms = pos-S;

smxn2 = sum(xms.^2,2);
smxn = sqrt(smxn2);
xn2 = sum(pos.^2,2);
xn = sqrt(xn2);

r = smxn + xn; % bistatic range 

ph = pos./xn;
xmsh = xms./smxn;
elipNorm = ph + xmsh; % normal to ellipsoid 
% velp = permute(vel,[2,1,3]);
% d = mm3d(elipNorm,velp); % bistatic doppler


% derivative of doppler with respect to position:
% term1 = mm3d(xms,velp).*xms./(smxn.*smxn2);
% term2 = mm3d(pos,velp).*pos./(xn.*xn2);
% dDdX = vel./smxn + vel./xn - term1 - term2;


if isempty(varargin) % normalize or not 
  y = r;
  jac = elipNorm;
else
  y = r./varargin{1};
  jac = elipNorm./varargin{1};
end
% y = perm(y,[2,1,3]);
