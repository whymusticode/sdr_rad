function xyzDet = disambigGrid(det)
% memory is what makes this concept prohibitive. So split the grid up



xx = gpuArray(-5e4:1e2:5e4-1);
yy = xx;
zz = gpuArray(0:2e2:20e3-1);
xyzDet = [];
for j = 1:10
  for k = 1:10
    for l = 1:10
      
      
      x = xx((1:100) + (j-1)*100);
      y = yy((1:100) + (k-1)*100);
      z = zz((1:10) + (l-1)*10);
      
      [X,Y,Z] = ndgrid(x,y,z);
      XYZ = [X(:),Y(:),Z(:)];
      R = sqrt(sum(XYZ.^2,2));
      matches = zeros(size(XYZ,1),1,'single','gpuArray');
      for i = 1:length(det)
        s = det(i).s;
        R2 = sqrt(sum((XYZ-s).^2,2));
        R2R = R2+R;
        idx = abs(R2R - det(i).rrd(:,1)') < 5e1;
        matches = matches + sum(idx,2);
      end
      idxDet = matches >= 4;
      
      xyzDet = [xyzDet;XYZ(idxDet,:)];
      
      
    end
  end
end
