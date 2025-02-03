

Nfast = 2^9;
Nslow = 2^11;
Ncpi = 2^5;
Nch = 2^2;
Nacc = 2^5;
rBin = 60;
vBin = 1;
aBin = 1;
Tcpi = .2;

r = (0:Nfast-1).*rBin;
v = (-Nslow/2:Nslow/2-1).*vBin;
a = (-Nacc/2:Nacc/2-1).*aBin;
t = (1:Ncpi).*Tcpi;
t = t - mean(t);
r2 = gS(r);
v2 = permute(gS(v),[2,1]);
a2 = permute(gS(a),[3,2,1]);
t2 = gS(t);

rva = Nfast*Nslow*Nacc;

data =  rand(Nfast*Nslow,Ncpi,Nch,'single','gpuArray');

nci(data,r2,v2,a2,t2,Nfast,Nslow,rBin,vBin,Ncpi)

% in = randn(512,2048,32,4,'single','gpuArray');
% ref = randn(512,1,'single','gpuArray');
% rd0 = @(in,ref)fft2(fft(in).*ref);
% tic 
% rd = rd0(ref,in);
% toc





