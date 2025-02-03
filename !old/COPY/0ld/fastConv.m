clear
m = 2e4;
n = 40;
r = 500; % number of actual bins that need to be computed

% x = gpuArray(randn(m,n));
% y = gpuArray(randn(m,n));
% x = randn(m,n);
% y = randn(m,n);
% algorithm 1
% z1 = gpuArray(zeros(2*r+1,n));
% z1 = zeros(2*r+1,n);

% tic
% for i = 1:n
%   z1(:,i) = xcorr(x(:,i),y(:,i),r);
% end
% z1 = z1(r+1:end,:);
% toc


% algorithm 2
compTime = 0;
allocTime = 0;
 x = gpuArray(randn(m,n));
  y = gpuArray(randn(m,n));
% x = randn(m,n);
%   y = randn(m,n);
for i = 1:1e0
  tic
  x = gpuArray(randn(m,n));
  y = gpuArray(randn(m,n));
%   x = randn(m,n);
%   y = randn(m,n);

  allocTime = allocTime + toc;
  tic
  for j = 1:1e3
  z2 = ifft(fft(x,2*m).*fft(flipud(y), 2*m));
  z2 = z2(m:m+r,:);
  end
  compTime = compTime+toc;
end
disp(compTime+allocTime)