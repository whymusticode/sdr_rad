clear
p = .9;
m = 8;
n = 6;
sig = .1;

x = randn(n,1);

H = linspace(-1,1,m)'.^(0:n-1);

yt = H*x;

R = eye(m)*sig^2;
W = R^-1;
P = (H'*W*H)^-1;
S = R - H*P*H';
% Si2 = S^-1;
Si = diag(1./diag(S));

res = [];
M = [];
for i = 1:1e3
  ym = yt + sig*randn(size(yt));
  xh = P*H'*W*ym ;
  
  yh = H*xh;
  
  res(:,i) = yh - ym;
  M2(i) = res(:,i)'*Si*res(:,i);
%   M3(i) = res(:,i)'*W*res(:,i);
%   tmp = pcg(S,res(:,i));
%   M4(i) = res(:,i)'*tmp;
end

thresh = chi2inv(p,m);

sum(M2 < thresh)/length(M2)
% sum(M3 < thresh)/length(M3)
% sum(M4 < thresh)/length(M4)


% cov(res')
% S


% chi2cdf(100,m)
% 
% 
% 
% 
% histogram(abs(M)-sqrt(m))
% disp(S)
% disp(cov(res'))

% sqrt(res'*R^-1*res)

% sqrt(res(:,1)'*S^-1*res(:,1))

% R - H*P*H' == cov(yh - y)

