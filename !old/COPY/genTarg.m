function out = genTarg(n)
txyz = [(rand(n,2)-.5)*1e5,5e3*ones(n,1)+5e3*ones(n,1)]; % target
vxyz = [randn(n,2)*100,zeros(n,1)]; % target velocities

out = [txyz';vxyz'];