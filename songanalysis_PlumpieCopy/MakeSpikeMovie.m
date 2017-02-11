function [] = MakeSpikeMovie(Spikes, Fs, OutputFileName)

FrameRate = 10;
Time = 0:1/Fs:length(Spikes)/Fs;
Time(end) = [];

NoofFrames = floor(length(Spikes)/(Fs/FrameRate));

for i = 1:NoofFrames,
    SI = 1;
    EI = i*(Fs/FrameRate);
    plot(Time(SI:EI), (Spikes(SI:EI)), 'k', 'LineWidth', 2);
    axis([0 Time(end) -5 5]);
    set(gca, 'Box', 'off');
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