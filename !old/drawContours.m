function [x,y,matches] = drawContours(det)

x = -100:.3:100;
y = x;
z = 0;

[X,Y,Z] = ndgrid(x,y,z);
XYZ = [X(:),Y(:),Z(:)];
R = sqrt(sum(XYZ.^2,2));
matches = zeros(size(XYZ,1),1,'single','gpuArray');
for i = 1:length(det)
  s = det(i).s/1e3;
  R2 = sqrt(sum((XYZ-s).^2,2));
  R2R = R2+R;
  idx = abs(R2R - det(i).rrd(:,1)'/1e3) < 1;
  matches = matches + sum(idx,2);
end

matches = reshape(matches,[length(y),length(x),length(z)])';



% idxDet = matches >= 4;
% 
% xyzDet = [xyzDet;XYZ(idxDet,:)];


