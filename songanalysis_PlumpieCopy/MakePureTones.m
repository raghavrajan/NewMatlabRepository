function [PureToneSeries] = MakePureTones(Freqs, Fs, Dur, InterToneDur)

t = (1:round(Dur*Fs))/Fs;
PureToneSeries = [zeros(1,(round(InterToneDur*Fs)))];
for i = 1:length(Freqs),
    PureTone = zeros(size(t));
    s = sin(2*pi*Freqs(i)*t);
    % include an on and off ramp of 1/10th the duration of the tone.
    OnRamp = cos(linspace(pi, 2*pi, Fs*Dur/10));
    OnRamp = OnRamp + 1;
    OnRamp = OnRamp/2;
    OffRamp = fliplr(OnRamp);

    r = [OnRamp(:)' ones(1, length(t) - length(OnRamp)*2) OffRamp(:)'];
    s = s.*r;
    PureTone = PureTone + s;
    PureTone = PureTone/(max(PureTone)*1.1);
    PureToneSeries = [PureToneSeries PureTone zeros(1,(round(InterToneDur*Fs)))];
end
audiowrite([num2str(Freqs(1)), '.to.', num2str(Freqs(end)), 'Hz.Dur.', num2str(Dur*1000), 'ms.InterToneDur.', num2str(InterToneDur*1000), '.PureToneSeries.wav'], PureToneSeries, Fs);
PlotSpectrogram(pwd, [num2str(Freqs(1)), '.to.', num2str(Freqs(end)), 'Hz.Dur.', num2str(Dur*1000), 'ms.InterToneDur.', num2str(InterToneDur*1000), '.PureToneSeries.wav'], 'wav', 'hot');