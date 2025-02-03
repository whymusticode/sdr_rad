clear
nTx = 5;
syms x y z real
exyz = sym('E', [nTx 3],'real');
rMeas = sym('r', [nTx 1],'real');
tg = [x,y,z];

r = sqrt(sum((exyz-tg).^2,2)) + sqrt(sum(tg.^2,2)); 

J = sum((r-rMeas).^2);

% syms rErr positive
% sig = ones(size(rdoa))*rErr;
% CRB = getCrb(rdoa,tg,[],sig);
% CRB = simplify(CRB);
% crb = matlabFunction(CRB);
% solu = solve(J==0,tg,'IgnoreAnalyticConstraints', true);

dJdx = simplify(diff(J,x));
dJdy = simplify(diff(J,y));
dJdz = simplify(diff(J,z));

sol = solve([dJdx,dJdy,dJdz]==[0,0,0],tg, ...
  'IgnoreAnalyticConstraints', true);

% simplify(S, 'Steps', 50)
% sol = solve(dJdx==0,x, 'IgnoreAnalyticConstraints', true);
