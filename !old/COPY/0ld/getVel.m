function vel = getVel(x,s,d)
 % least squares bistatic cartesian velocity given position 
% INPUT:
% x is [1x3xN] position of target
% s is [Mx3xN] position of emitters 
% d is [Mx1xN] doppler measured by each emitter 
%
% OUTPUT:
% vel is [1x3xN] velocity of target

xms = x - s;

normx = sqrt(sum(x.^2,2));
normxms = sqrt(sum(xms.^2,2));

xmsh = xms./normxms;
xh = x./normx;

H = xmsh + xh; % [Mx3xN]


Hp = permute(H,[2,1,3]);

H2 = mm3d(Hp,H);

H2i = pagefun(@inv,H2);

vel = mm3d(mm3d(H2i,Hp),d);

vel = permute(vel,[2,1,3]);