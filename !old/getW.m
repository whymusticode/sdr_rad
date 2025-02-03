function eopt = getW(Ntaps,tx_signal,rx_signal) %RLSest,LMSest,
N = length(tx_signal);
% R = 0;
% P = 0;
% nnn = 0;
% xvec = zeros(Ntaps,1);
% for nn=1:N
%   nnn = nnn+1;
%   xvec = circshift(xvec,1);
%   xvec(1) = tx_signal(nn);
%   R = R+xvec*xvec';
%   P = P + conj(d(nn))*xvec;
%   
%   
% end
% R = R/nnn;
% P = P/nnn;
% 
% wopt = pinv(R)*P;       %optimal filter weights to obtain MMSE
% 
% xvec = zeros(Ntaps,1);
% for nn=1:N
%    xvec = circshift(xvec,1);
%   xvec(1) = tx_signal(nn);
%   est(nn) = wopt'*xvec;
% %   eopt(nn) = d(nn) - wopt'*xvec;
% end



% rx_signal = gather(rx_signal);
% tx_signal = gather(tx_signal);


d = rx_signal;
%% Compute MMSE and optimal filter weights
n = length(tx_signal);
Hidx = (1-Ntaps:n-Ntaps)'+(1:Ntaps);
idx = Hidx<1 | Hidx > n;
Hidx(idx) = 1;
H = tx_signal(Hidx);
H(idx) = 0;
wopt = (H'*H)^-1*H'*d;

eopt = d-H*wopt;

plot(1:n,real(d),1:n,real(H*wopt))

% imagesc(abs(R))

% R = 0;
% P = 0;
% nnn = 0;
% xvec = zeros(Ntaps,1);
% for nn=1:N
%   nnn = nnn+1;
%   xvec = circshift(xvec,1);
%   xvec(1) = tx_signal(nn);
%   R = R+xvec*xvec';
%   P = P + conj(d(nn))*xvec;
% end
% R = R/nnn;
% P = P/nnn;
% 
% wopt = pinv(R)*P;       %optimal filter weights to obtain MMSE
% return
% yopt = filter(conj(wopt),1,tx_signal);
% eopt = d-yopt;

% figure;
% hold on;
% plot(db(d));
% plot(db(eopt));
% return


%% RLS and LMS adaptive filters
% eopt = zeros(N,1);
% xvec = zeros(Ntaps,1);
% 
% %LMS parameters
% muLMS = .05/trace(R);               
% %this can be tuned to control LMS convergence and stability....
% wLMS = zeros(Ntaps,1);
% eLMS = zeros(N,1);
% 
% %RLS parameters
% deltaRLS = 1e-3;                    
% %this can be tuned to control RLS convergence and stability....
% lambdaRLS = .99;                    
% %this can be tuned to control RLS convergence and stability....
% wRLS = zeros(Ntaps,1);
% eRLS = zeros(N,1);
% eRLS_mat = zeros(N,1);
% InvpsiRLS = eye(Ntaps,Ntaps)/deltaRLS;
% 
% %built-in MATLAB RLS filter for comparison
% h_rls = dsp.RLSFilter(Ntaps, 'ForgettingFactor', lambdaRLS, ...
%   'InitialInverseCovariance', InvpsiRLS);


% for nn=1:N
%   xvec = circshift(xvec,1);
%   xvec(1) = tx_signal(nn);
%   
%   yLMS = wLMS'*xvec;               %y(n) is output from adaptive filter
%   eLMS(nn) = d(nn) - yLMS;        %e(n) is error signal
%   wLMS = wLMS + 2*muLMS*conj(eLMS(nn))*xvec;       
%   %w is adpative weight vector computed by LMS algorithm
%   
%   
%   g  = InvpsiRLS*xvec;
%   k = g/(lambdaRLS+xvec'*g);
%   yRLS = wRLS'*xvec;               %y(n) is output from adaptive filter
%   eRLS(nn) = d(nn) - yRLS;        %e(n) is error signal
%   wRLS = wRLS + k*conj(eRLS(nn));       
%   %w is adpative weight vector computed by RLS algorithm
%   InvpsiRLS = (eye(Ntaps)-k*xvec')*InvpsiRLS/lambdaRLS;
%   
%   [~,eRLS_mat(nn)] = step(h_rls, tx_signal(nn), d(nn));
%   
%   eopt(nn,1) = d(nn) - wopt'*xvec;   
%   %keep track of resulting error from MMSE fitler as well
%   
%   RLSest(nn) = yRLS;
%   LMSest(nn) = yLMS;
  
%   opt(nn,1) = d(nn) - wopt'*xvec;
  
  
% end


% plot(1:nn,real(opt),1:nn,real(d))



