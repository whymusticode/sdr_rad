function idx = uniquify(lam,direc)
% returns an indexing scheme that returns  indicies of the strongest towers
% at each wavelength 
% 
% 
% lam is a vector of wavelengths
% direc is a vector of power recieved 


[DIR,sidx] = sort(direc);
id2 = ~isnan(DIR);
% DIR = DIR(id2);
sidx = sidx(id2);
la = lam(sidx);


[~,lastDX,~] = unique(la,'last');

idx = flipud(sidx(lastDX));

% step 2, look for any signals drowned out in sidelobes
d = direc(idx);

f = 3e2./lam(idx);


sidelobeIDX = ~any(abs(f-f') <= .4 & d' - d > 20,2);
 
idx = idx(sidelobeIDX);

















