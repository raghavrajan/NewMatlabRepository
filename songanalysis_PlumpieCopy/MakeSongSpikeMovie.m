function [] = MakeSongSpikeMovie(Song, Spikes, Fs, OutputFileName)

MaxSpectVal = CalculateMaxSpectVal(Song, Fs);

FrameRate = 25;
Time = 0:1/Fs:length(Song)/Fs;
Time(end) = [];

NoofFrames = floor(length(Song)/(Fs/FrameRate));

VideoObj = VideoWriter(OutputFileName);
VideoObj.FrameRate = FrameRate;

open(VideoObj);

SpikeMultiplier = 20;

for i = 1:NoofFrames,
    SI = 1;
    EI = i*(Fs/FrameRate);
    PlotSpectrogram_SongVar_Scaled(Song(SI:EI), (1:1:length(Song(SI:EI)))/Fs, Fs, MaxSpectVal)
    %set(gcf, 'Position', [680 558 800 450]);
    hold on;
    plot(Time(SI:EI), (Spikes(SI:EI)*SpikeMultiplier)+ 10000 + 1000 + abs(min(Spikes*SpikeMultiplier)), 'k', 'LineWidth', 0.5);
    axis([0 Time(end) 300 1.01*(max(Spikes*SpikeMultiplier) + 10000 + 1000 + abs(min(Spikes*SpikeMultiplier)))]);
    set(gcf, 'Color', 'w');
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    set(gca, 'Visible', 'off');
    F(i) = getframe(gcf);
    writeVideo(VideoObj, F(i));
    close all;
end
close(VideoObj);
movie(figure, F, 1, 10);

% movie2avi(F, OutputFileName, 'compression', 'None', 'fps', FrameRate);