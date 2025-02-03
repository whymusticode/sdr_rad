function   [x0,y0] = subP(x,y,X,Y,im)



x0 = [];
y0 = [];
for i = 1:length(X)
  if abs(y(Y(i))) < 5
    continue
  end
%   try
    ydx = Y(i)+(-1:1);
    xdx = X(i)+(-1:1);
    in = gather(im(ydx,xdx));
%     x0(end+1) = x(xdx(2));
%     y0(end+1) = y(xdx(2));
    [x0(end+1),y0(end+1)] = getSubPixel3(x(xdx),y(ydx),in);
%   end
  
  
  
end









