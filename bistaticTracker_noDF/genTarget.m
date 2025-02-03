function x = genTarget(p)
% inits a new target

side = randi(4);
th = (side+1)*pi/2 - pi/4 + pi/2*rand();

switch side
  case 1
    x = [p.lim,(rand()-.5)*p.lim];
  case 2
    x = [(rand()-.5)*p.lim,p.lim];
  case 3
    x = [-p.lim,(rand()-.5)*p.lim];
  case 4
    x = [(rand()-.5)*p.lim,-p.lim];
end

x = [x,5e3+1e4*rand(),(rotMat(-th)*[300;0])',0];
