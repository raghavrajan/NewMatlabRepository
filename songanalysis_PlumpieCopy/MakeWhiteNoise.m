function [Output] = MakeWhiteNoise(TotalDur, PlotToggle)

Fs = 44100;
Output = rand(round(44100*TotalDur),1);
audiowrite('WhiteNoise.wav', Output, Fs);

if (strfind(PlotToggle, 'on'))
    PlotSpectrogram_SongVar(Output, Fs);
end