function y = Hbistat(pos,vel,S,varargin)
% strange format

% output is [emitters,range then doppler,targets]

% pos = permute(posVel(:,1:3),[3,2,1]);
% vel = permute(posVel(:,4:6),[3,2,1]);

% pos 1x3xN
% vel ''
% s   Mx3xN

xms = pos-S;

smxn = sqrt(sum(xms.^2,2));
xn = sqrt(sum(pos.^2,2));

t = smxn + xn;

ph = pos./xn;
xmsh = xms./smxn;
tdi = mm3d((ph+xmsh),permute(vel,[2,1,3]));

if isempty(varargin)
  y = [t,tdi];% return non- normalized values
else
  y = [t./varargin{1},tdi./varargin{2}];
  % normalize this for iterative min solver
end

