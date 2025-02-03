clear
nE = 3;
syms x y xd yd zd real
syms sigd sigr z positive
exy = sym('E', [nE 2],'real');
ez = sym('E', [nE 1],'positive');
exyz = [exy,ez];
rMeas = sym('rr', [nE 1],'positive');
dMeas = sym('dd', [nE 1],'positive');

tg = [x,y,z];
d = [xd,yd,zd];

xms = exyz-tg;
normXms = sqrt(sum(xms.^2,2));
normr = sqrt(sum(tg.^2,2));
r = normXms + normr; 
rd = (xms./normXms + tg./normr)*d';
% taylor(rd,tg)
% assume(r,'positive')
% eqnD = [simplify(rd == dMeas);r==rMeas];
% determinedSolution = solve(eqnD,d,'IgnoreAnalyticConstraints', true);
% 
% detSol.xd = simplify(determinedSolution.xd);
% detSol.yd = simplify(determinedSolution.yd);
% detSol.zd = simplify(determinedSolution.zd);


eqnRD = simplify([rd == dMeas;r==rMeas]);
determinedSolution = solve(eqnRD,[tg,d],'IgnoreAnalyticConstraints', true);

return


J = sum([(t-rMeas).^2/sigr^2;(td-dMeas).^2/sigd^2]);

J = sum((td-dMeas).^2);
dJdxd = simplify(diff(J,xd));
dJdyd = simplify(diff(J,yd));
dJdzd = simplify(diff(J,zd));
sol = solve([dJdxd,dJdyd,dJdzd]==zeros(3,1),d, ...
  'IgnoreAnalyticConstraints', true);



assume(J,'positive')
dJdx = simplify(diff(J,x));
dJdy = simplify(diff(J,y));
dJdz = simplify(diff(J,z));
dJdxd = simplify(diff(J,xd));
dJdyd = simplify(diff(J,yd));
dJdzd = simplify(diff(J,zd));
sol = solve([dJdx,dJdy,dJdz,dJdxd,dJdyd,dJdzd]==zeros(6,1),[r,rd], ...
  'IgnoreAnalyticConstraints', true);

% simplify(S, 'Steps', 50)
% sol = solve(dJdx==0,x, 'IgnoreAnalyticConstraints', true);


%% 1d case has solutions
% syms t x s
% solve(t==sqrt((x-s)^2)+sqrt(x^2))
%  s/2 - t/2
%  s/2 + t/2

%% 2d case has no solution
% syms x y s11 s12 s21 s22 real
% syms t1 t2 positive
% 
% sol2d = solve([t1==sqrt((x-s11)^2+(y-s12)^2)+sqrt(x^2+y^2),...
%   t2==sqrt((x-s21)^2+(y-s22)^2)+sqrt(x^2+y^2)],x,...
%   'IgnoreAnalyticConstraints', true);


%% this is look at the cost function for the 2d case 
% clear
% % syms x y z real
% % exyz = sym('E', [5 2],'real');
% % rdoaMeas = sym('r', [5 1],'real');
% tg = [-1.5,0];
% 
% exyz = randn(5,2);
% rd = sqrt(sum((exyz-tg).^2,2)) + sqrt(sum(tg.^2,2)); 
% 
% x = linspace(-2,2,100);
% y = x;
% [X,Y] = meshgrid(x,y);
% 
% 
% rdMeas = rd + 1e-6*randn(5,1);
% test = permute([X(:),Y(:)],[3,2,1]);
% 
% rdtest = sqrt(sum((exyz-test).^2,2)) + sqrt(sum(test.^2,2)); 
% J = sum((rdtest - rdMeas).^2);
% J = reshape(J,size(X));
% figure(34)
% hold off
% imagesc(x,y,J)
% hold on
% plot(tg(1),tg(2),'wo')
% miniMask = J<circshift(J,[1,0])&J<circshift(J,[-1,0])&...
%   J<circshift(J,[0,1])&J<circshift(J,[0,-1])& ...
%   J<circshift(J,[1,1])&J<circshift(J,[-1,-1])&...
%   J<circshift(J,[-1,1])&J<circshift(J,[1,-1]);
% plot(X(miniMask),Y(miniMask),'w+')
% plot(exyz(:,1),exyz(:,2),'wd')
% test = permute(gpuArray([X(:),Y(:)]),[1,3,2]);













