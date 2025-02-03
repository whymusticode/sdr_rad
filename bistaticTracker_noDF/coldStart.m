function y = coldStart(S)
% [y,idx] = TDOA2XYZ(S)
% [Xest, Pest] = TDOA2XYZ(S,W)
% this code finds the solution to S(:,4) = dis(X,S(:,1:3)) + b
% where S = [x1,y1,z1,t1;x2,y2 ... ] and dis() is distance function
% the x y and z values are position values of sensor locations t1, t2,
% etc... are position differences relative to one of the sensors, for
% example, if you compute the time difference of arival between 3 and 1,
% call that td31, t3 = c*td31 (c is the speed of light (or more generally
% the thing measured)) t1 = 0, t2 = c*td21, t3 = c*td31, etc... S NEEDS TO
% HAVE >= 4 rows this is based on "An Algebraic Solution of the GPS
% Equations" by Stephen Bancroft of King Radio

mk = @(a,b) sum(a.*b.*repmat([1,1,1,-1],[size(a,1),1,size(S,3)]),2);
% ^ minkowski functional

r = mk(S,S)/2;
Sp = permute(S,[2,1,3]);
B = ml3d(mm3d(Sp,S),Sp); 
% times rms^2 for true covar

u = sum(B,2);
v = mm3d(B,r);
up = permute(u,[2,1,3]);
vp = permute(v,[2,1,3]);

E = mk(up,up);
F = mk(up,vp)-1;
G = mk(vp,vp);

f2eg = F.^2-E.*G; % most of these are < 0 indicating no solution 
idx = f2eg < 0; 
% about 95% pd and 5% pfa for using this as a metric to remove ghosts
f2eg(idx) = 0;

sfeg = sqrt(f2eg);
fidx = -F;
lam = [fidx+sfeg,fidx-sfeg]./E;

yPre = mm3d(u,lam)+repmat(v,[1,2,1]); 


% just some minor formatting  only return a solution z > 0
y = yPre(:,1,:);
predx = y(3,:,:) < 0;
y(:,:,predx) = yPre(:,2,predx);
y(4,:,:) = - y(4,:,:);

% only return a solution z > 0

% P = covar(:,:,idx);


% ^ two potentially correct estimates per input, note both solutions will
% have similar accuracy

% idx(linspace(1,5^7,5))
% tmp = y(4,:,:);
% histogram(y(4,:,:))

% this part removes one incorrect estimate, by plugging back in
% err = [0,0];
% for i = 1:2
%   err(i) = abs(dis(y(1:dim-1,i),S(1,1:dim-1)')-y(dim,i) - S(1,dim));
% end
% [~,loc] = min(err);
% Xest = y(1:dim-1,loc)';
% err = err(loc);

% varargout = {};
% if ~getCovar
%   return
% end
% % covar = (S(:,1:3)'*W*S(:,1:3))^-1;
% varargout = {covar};



% y = simplify(y);




%% this allows the operation to be computed symbolically 
%  (place at beginning of code)
% clear
% nSen = 4;
% sMat = '[';
% for i = 1:nSen
%   strI = num2str(i);
%   thing = ['x',strI,' y',strI,' z',strI,' t',strI];
%   temp = evalc(['syms ',thing,' real']);
%   sMat = [sMat,thing,';'];
% end
% sMat = [sMat(1:end-1),']'];
% temp = evalc(['S = ',sMat]);
