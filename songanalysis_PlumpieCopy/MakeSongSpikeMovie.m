function [] = MakeSongSpikeMovie(Song, Spikes, Fs, OutputFileName)

MaxSpectVal = CalculateMaxSpectVal(Song, Fs);

FrameRate = 10;
Time = 0:1/Fs:length(Song)/Fs;
Time(end) = [];

NoofFrames = floor(length(Song)/(Fs/FrameRate));

for i = 1:NoofFrames,
    SI = 1;
    EI = i*(Fs/FrameRate);
    PlotSpectrogram_SongVar_Scaled(Song(SI:EI), Fs, MaxSpectVal)
    hold on;
    plot(Time(SI:EI), (Spikes(SI:EI)*1500)+ 20000, 'k', 'LineWidth', 2);
    axis([0 Time(end) 300 32000]);
    set(gcf, 'Color', 'w');
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    set(gca, 'Visible', 'off');
    F(i) = getframe(gcf);
    F(i) = getframe;
    c;
end
movie(figure, F, 1, 10);
movie2avi(F, OutputFileName, 'compression', 'None', 'fps', FrameRate);