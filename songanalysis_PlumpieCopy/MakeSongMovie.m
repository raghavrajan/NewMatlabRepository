function [] = MakeSongMovie(Song, Fs, OutputFileName)

MaxSpectVal = CalculateMaxSpectVal(Song, Fs);

FrameRate = 10;
Time = 0:1/Fs:length(Song)/Fs;
Time(end) = [];

NoofFrames = floor(length(Song)/(Fs/FrameRate));

for i = 1:NoofFrames,
    SI = 1;
    EI = i*(Fs/FrameRate);
    PlotSpectrogram_SongVar_Scaled(Song(SI:EI), Fs, MaxSpectVal)
    axis([0 Time(end) 300 8000]);
    set(gcf, 'Color', 'w');
    set(gca, 'XTick', [0:0.2:Time(end)]);
    set(gca, 'FontSize', 20, 'FontName', 'Arial');
    set(gca, 'YTick', []);
    set(gca, 'YColor', 'w');
    set(gca, 'Visible', 'off');
    F(i) = getframe(gcf);
    F(i) = getframe;
    c;
end
movie(figure, F, 1, 10);
movie2avi(F, OutputFileName, 'compression', 'None', 'fps', FrameRate);