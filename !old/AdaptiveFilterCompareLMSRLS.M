% Adaptive Filter Comparison

close all;
clear all;

Ntaps = 16; %number of filter taps

%N = N;
N = 3000;

tx_signal = randn(N,1) + 1i*randn(N,1);

%% Shape Tx spectrum
%set to 0 for unfiltered waveform (white noise, full digital spectrum)
%set to 1 for band limited waveform (pink noise)

%for simple LMS, white noise converges well. Bandlimited waveform has
%convergence problems due to ill-conditioned performance surface
if(1)
    FilterOrder = 256;
    FilterF1 = .3;
    FilterF2 = .7;
    d = fdesign.bandpass('N,F3dB1,F3dB2',FilterOrder,FilterF1,FilterF2);
    f_BPF = design(d); %,'equiripple');
    tx_signal = filter(f_BPF,tx_signal(:)); %Filtered Tx signal
end;
tx_signal = tx_signal/max(tx_signal); %normalize

%% set up receive signal
% rx signal is filtered through random channel
channel = randn(Ntaps,1) + 1i*randn(Ntaps,1);
rx_signal = filter(channel,1,tx_signal); %Received signal
rx_signal = rx_signal/max(rx_signal); %normalize

rx_signal = rx_signal(1:N);

d = rx_signal;


addnoise = 1;
if(addnoise)
    txnoisedB = -75;
    dnoisedB = -75;
    
    d = d+randn(size(d))*10^(dnoisedB/20)*sqrt(var(d));
    tx_signal = tx_signal+randn(size(tx_signal))*10^(txnoisedB/20)*sqrt(var(tx_signal));
end

%% Compute MMSE and optimal filter weights
R = 0;
P = 0;
nnn = 0;
xvec = zeros(Ntaps,1);
for nn=1:N
    nnn = nnn+1;
    xvec = circshift(xvec,1);
    xvec(1) = tx_signal(nn);
    R = R+xvec*xvec';
    P = P + conj(d(nn))*xvec;
end
R = R/nnn;
P = P/nnn;

wopt = pinv(R)*P;       %optimal filter weights to obtain MMSE

% yopt = filter(conj(wopt),1,tx_signal);
% eopt = d-yopt;

% figure;
% hold on;
% plot(db(d));
% plot(db(eopt));
% return


%% RLS and LMS adaptive filters
eopt = zeros(N,1);
xvec = zeros(Ntaps,1);

%LMS parameters
muLMS = .05/trace(R);               %this can be tuned to control LMS convergence and stability....
wLMS = zeros(Ntaps,1);
eLMS = zeros(N,1);

%RLS parameters
deltaRLS = 1e-3;                    %this can be tuned to control RLS convergence and stability....
lambdaRLS = .99;                    %this can be tuned to control RLS convergence and stability....
wRLS = zeros(Ntaps,1);
eRLS = zeros(N,1);
eRLS_mat = zeros(N,1);
InvpsiRLS = eye(Ntaps,Ntaps)/deltaRLS;

%built-in MATLAB RLS filter for comparison
h_rls = dsp.RLSFilter(Ntaps, 'ForgettingFactor', lambdaRLS, 'InitialInverseCovariance', InvpsiRLS);


for nn=1:N
    xvec = circshift(xvec,1);
    xvec(1) = tx_signal(nn);

    yLMS = wLMS'*xvec;               %y(n) is output from adaptive filter
    eLMS(nn) = d(nn) - yLMS;        %e(n) is error signal
    wLMS = wLMS + 2*muLMS*conj(eLMS(nn))*xvec;       %w is adpative weight vector computed by LMS algorithm

    
    g  = InvpsiRLS*xvec;
    k = g/(lambdaRLS+xvec'*g);
    yRLS = wRLS'*xvec;               %y(n) is output from adaptive filter
    eRLS(nn) = d(nn) - yRLS;        %e(n) is error signal
    wRLS = wRLS + k*conj(eRLS(nn));       %w is adpative weight vector computed by RLS algorithm
    InvpsiRLS = (eye(Ntaps)-k*xvec')*InvpsiRLS/lambdaRLS;
    
    [~,eRLS_mat(nn)] = step(h_rls, tx_signal(nn), d(nn));
    
    eopt(nn) = d(nn) - wopt'*xvec;   %keep track of resulting error from MMSE fitler as well

    
end


% figure;
% hold on;
% plot(real(wLMS)/max(abs(wLMS)));
% plot(real(wRLS)/max(abs(wRLS)));
% plot(real(wopt)/max(abs(wopt)));
% 
% figure;
% hold on;
% plot(imag(wLMS)/max(abs(wLMS)));
% plot(imag(wRLS)/max(abs(wRLS)));
% plot(imag(wopt)/max(abs(wopt)));


% 
%  hftf1 = dsp.FastTransversalFilter(Ntaps,'ForgettingFactor',.98);
%  [~,eFTRLS] = step(hftf1,tx_signal,d);

figure;
hold on;
plot(absdB(d));
plot(absdB(eLMS));
plot(absdB(eRLS));
plot(absdB(eRLS_mat));
plot(absdB(eopt));
axis([1 N, -100, 0])
ylabel 'Error (dB)';
xlabel 'Samples';
legend('Received', 'LMS', 'RLS', 'MATLAB RLS', 'MMSE');
title 'Error Convergence'
%boldify;


figure;
hold on;
Nfft = 1024;
f = linspace(0,1,Nfft);
plot(f, absdB(fft(tx_signal(end-Nfft+1:end))));
plot(f, absdB(fft(d(end-Nfft+1:end))));
plot(f, absdB(fft(eLMS(end-Nfft+1:end))));
plot(f, absdB(fft(eRLS(end-Nfft+1:end))));
plot(f, absdB(fft(eRLS_mat(end-Nfft+1:end))));
plot(f, absdB(fft(eopt(end-Nfft+1:end))));
%boldify;
axis([0 .5, -70, 40])
ylabel 'dB';
xlabel 'Normalized Freq (Fs)';
title 'Canceller Spectrum';
legend ('Transmitted', 'Received', 'LMS', 'RLS', 'RLS MATLAB', 'MMSE');


