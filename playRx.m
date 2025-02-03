function playRx(rcv, fs)
% plays recorded data

% n = length(rcv);
% steps = floor(n/7e3);
% rcv = rcv(1:steps*7e3);

freqDeviation = 75000;
filterTimeConst = 7.5e-5;

audioSampleRate = 48000;

fmBroadcastDemod = comm.FMBroadcastDemodulator(...
  'SampleRate', fs,'FrequencyDeviation', freqDeviation, ...
  'FilterTimeConstant', filterTimeConst,...
  'AudioSampleRate', audioSampleRate, 'Stereo', true,'RBDS',true);

player = audioDeviceWriter('SampleRate',audioSampleRate);

for i = 1:size(rcv,2)
[audioSig,rcvdm] = fmBroadcastDemod(rcv(:,i)); %#ok<ASGLU>

player(audioSig);
end


%% end used code


% radioTime = 0;
% frontEndFrameTime = samplesPerFrame/fs;

%%
% Create an RTL-SDR receiver System object.
% fmReceiver = comm.SDRRTLReceiver('OutputDataType','double',...
%     'SampleRate',fs,'SamplesPerFrame', samplesPerFrame,...
%     'CenterFrequency',fc);
%%
% Create an FM broadcast receiver System object.


