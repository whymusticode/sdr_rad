function playRx(rcv, fs)
% plays recorded data



% fc = station*1e6; % FM radio station
% samplesPerFrame = 3840;
% fs = 240000;

freqDeviation = 75000;
filterTimeConst = 7.5e-5;

audioSampleRate = 48000;

% radioTime = 0;
% frontEndFrameTime = samplesPerFrame/fs;

%%
% Create an RTL-SDR receiver System object.
% fmReceiver = comm.SDRRTLReceiver('OutputDataType','double',...
%     'SampleRate',fs,'SamplesPerFrame', samplesPerFrame,...
%     'CenterFrequency',fc);
%%
% Create an FM broadcast receiver System object.
fmBroadcastDemod = comm.FMBroadcastDemodulator(...
  'SampleRate', fs,...
  'FrequencyDeviation', freqDeviation, 'FilterTimeConstant', filterTimeConst,...
  'AudioSampleRate', audioSampleRate, 'Stereo', true,'RBDS',true);
%%
% Create an audio device writer System Object(TM).
player = audioDeviceWriter('SampleRate',audioSampleRate);

%%
% Loop until timer expires, receiving and playing the specified FM station.
% while radioTime < duration
%     %%
%     % Receive baseband samples from the specified FM station.
%     [rcv,~,lost,~] = fmReceiver();
%%
% Demodulate FM broadcast signals and play the decoded audio.
[audioSig,rcvdm] = fmBroadcastDemod(rcv);

%     z = fmdemod(rcv,centerFreq,rfSampleRate,freqDeviation);
%     plot(rcvdm)
%     plot(real(rcv))
player(audioSig);
%%
% Update radio time. If there were lost samples, add those too.
%     radioTime = radioTime + frontEndFrameTime + double(lost)/rfSampleRate;
% end
%
% release(fmReceiver)
% release(fmBroadcastDemod)
% release(player)