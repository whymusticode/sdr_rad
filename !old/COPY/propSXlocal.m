function SX = propSXlocal(SX,ts,Q)
SX(1:3,:) = ts*SX(4:6,:)+SX(1:3,:);
SX = SX + mvnrnd(zeros(size(SX,1),1),Q,size(SX,2))';
end