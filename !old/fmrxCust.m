function fmrxCust(duration, station)
% This function builds an FM stereo receiver using MATLAB(R) and Communications 
% System Toolbox(TM), and enables you to receive signals in real time using the 
% RTL-SDR USB radio hardware and the Communications System Toolbox Support Package 
% for RTL-SDR Radio. Function inputs specify the duration in seconds to play the FM 
% signal and the FM station to tune the radio to.

% Initialize run-time parameters.
centerFreq = station*1e6; % FM radio station
samplesPerFrame = 3840;
rfSampleRate = 240000;

freqDeviation = 75000;
filterTimeConst = 7.5e-5;

audioSampleRate = 48000;

radioTime = 0;
frontEndFrameTime = samplesPerFrame/rfSampleRate;

%%
% Create an RTL-SDR receiver System object.
fmReceiver = comm.SDRRTLReceiver('OutputDataType','double',...
    'SampleRate',rfSampleRate,'SamplesPerFrame', samplesPerFrame,...
    'CenterFrequency',centerFreq);
%%
% Create an FM broadcast receiver System object.
fmBroadcastDemod = comm.FMBroadcastDemodulator(...
    'SampleRate', rfSampleRate,...
    'FrequencyDeviation', freqDeviation, 'FilterTimeConstant', filterTimeConst,...
    'AudioSampleRate', audioSampleRate, 'Stereo', true,'RBDS',true);
%%
% Create an audio device writer System Object(TM).
player = audioDeviceWriter('SampleRate',audioSampleRate);

%%
% Loop until timer expires, receiving and playing the specified FM station.
while radioTime < duration
    %%
    % Receive baseband samples from the specified FM station.
    [rcv,~,lost,~] = fmReceiver();
    %%
    % Demodulate FM broadcast signals and play the decoded audio.
    [audioSig,rcvdm] = fmBroadcastDemod(rcv);
    
%     z = fmdemod(rcv,centerFreq,rfSampleRate,freqDeviation);
%     plot(rcvdm)
%     plot(real(rcv))
    player(audioSig);
    %%
    % Update radio time. If there were lost samples, add those too.
    radioTime = radioTime + frontEndFrameTime + double(lost)/rfSampleRate;
end

release(fmReceiver)
release(fmBroadcastDemod)
release(player)