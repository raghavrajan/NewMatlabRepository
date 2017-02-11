function [PureTone] = MakePureTone(Freq, Fs, Dur, TotalDur, Amplitude)

t = (1:round(Dur*Fs))/Fs;
PureTone = zeros(size(t));
s = sin(2*pi*Freq*t);
% include an on and off ramp of 1/10th the duration of the tone.
OnRamp = cos(linspace(pi, 2*pi, Fs*Dur/10));
OnRamp = OnRamp + 1;
OnRamp = OnRamp/2;
OffRamp = fliplr(OnRamp);

%r = sin(linspace(0,pi/2,round(length(t)/10)));
r = [OnRamp(:)' ones(1, length(t) - length(OnRamp)*2) OffRamp(:)'];
s = s.*r;
PureTone = PureTone + s;

PureTone = [zeros(1,(round(TotalDur*Fs) - length(t))/2) PureTone zeros(1,(round(TotalDur*Fs) - length(t))/2)];
PureTone = PureTone + rand(size(PureTone))/1000;
PureTone = PureTone/(max(PureTone)*1.1);
audiowrite([num2str(Freq), 'Hz.', num2str(Dur*1000), 'ms.', num2str(Amplitude), '.PureTone.wav'], PureTone, Fs);
PlotSpectrogram(pwd, [num2str(Freq), 'Hz.', num2str(Dur*1000), 'ms.', num2str(Amplitude), '.PureTone.wav'], 'wav', 'hot');