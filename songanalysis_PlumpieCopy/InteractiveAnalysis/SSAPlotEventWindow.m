function [] = SSAPlotEventWindow

figure();
set(gcf, 'Color', 'w');
for i = 1:length(FileInfo.WSpikeTrain),
    line([[FileInfo.WSpikeTrain{i}]'; [FileInfo.WSpikeTrain{i}]'], [(ones(size([FileInfo.WSpikeTrain{i}]))'*i + 0.25); (ones(size([FileInfo.WSpikeTrain{i}]))'*i - 0.25)], 'Color', 'k', 'LineWidth', 0.5);
end
axis tight;
hold on;
temp = axis;

for i = 1:length(EventParameters.EventTimes),
    plot([(EventParameters.EventTimes(i) - EventWindowWidth) (EventParameters.EventTimes(i) - EventWindowWidth)], [temp(3) temp(4)], 'b');
    plot([(EventParameters.EventTimes(i) + EventWindowWidth) (EventParameters.EventTimes(i) + EventWindowWidth)], [temp(3) temp(4)], 'b');
end
axis([-PreSongStartDuration MedianMotif.Length temp(3) temp(4)]);