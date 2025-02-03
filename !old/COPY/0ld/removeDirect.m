function retFull = removeDirect(nTaps,rx1,rx2) %RLSest,LMSest,
% Ntaps = 50;
n = length(rx1);

winLen = 1e2;



Hidx = (1-nTaps:n-nTaps)'+(1:nTaps) + round(nTaps/2);
idx = Hidx<1 | Hidx > n;
Hidx(idx) = 1;
H = rx1(Hidx);
H(idx) = 0;
% tic
retFull = zeros(size(rx1));
for i = 1:floor(n/winLen)
  subIdx = (1:winLen) +(i-1)*winLen;
  wopt = (H(subIdx,:)'*H(subIdx,:))^-1*H(subIdx,:)'*rx2(subIdx,:);
  ret = rx2(subIdx,:) - H(subIdx,:)*wopt;
  retFull(subIdx) = ret;
  % disp(10*log10(real(ret'*ret)/real(rx2(subIdx,:)'*rx2(subIdx,:))))
end

% plot(real(retFull))

%% RLS
% tic
% y = gather(rx2);
% H2 = gather(H);
% l = 1.04;
% X = zeros(Ntaps,1);
% I = eye(Ntaps);
% P = I;
% res = zeros(n,1);
% dir = zeros(n,1);
%
% for i = 1:n
%   h = H2(i,:);
%
%   P = P*l;
%
%   S = h*P*h' + 1;
%   dir(i) = h*X;
%   res(i) = y(i) - dir(i);
%
%   K = P*h'/S;
%   X = X + K*res(i);
%   P = (I-K*h)*P;
% end
%
% retu = res(101:5e4+100);
% direct = dir(101:5e4+100);
%
%
% plot(1:length(dir),real(dir),1:length(res),real(res))



% toc
% figure(2)
% plot(real(res))
% ylim([-1,1])
% grid on



% check indexes at the end
% H(1,:)
% H(end,:)



% retFull(1:5e4)'*retFull(1:5e4)
% res(1:5e4)'*res(1:5e4)

% plot(subIdx,real(rx2(subIdx,:)),subIdx,real(H(subIdx,:)*wopt))



% wopt = (H'*H)^-1*H'*rx2;

% eopt = rx2-H*wopt;

% real(eopt'*eopt)/real(rx2'*rx2)


% plot(1:n,real(rx2),1:n,real(H*wopt))
%
%
% plot(1:n,real(rx2),1:n,real(eopt))
%
% plot(real(wopt))




% [val,lag] = xcorr(rx1,rx2);
%   plot(lag,abs(val))