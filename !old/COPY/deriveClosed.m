



S = sym('S',[2,2]);
x = sym('X',[1,4]);
t = sym('t',[2,1]);
d = sym('d',[2,1]);
pos = x(1:2);
vel = x(3:4);

xms = pos-S;

smxn = sqrt(sum(xms.^2,2));
xn = sqrt(sum(x(1:3).^2,2));

ph = pos./xn;
xmsh = xms./smxn;


% d == (ph+xmsh)*vel';

sol = solve(t == smxn + xn,pos);

th(1) = simplify([sol.X1(1),sol.X2(1)]);
th(2) = simplify([sol.X1(2),sol.X2(2)]);





