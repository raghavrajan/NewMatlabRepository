function [] = SSAPlotExampleSpectST(DirectoryName, FileName, FileType, SpecAxis, SpikeTrainAxis, SpikeTimes, TimeRange, FileDur, PreTime, PostTime)

cla(SpecAxis);
axes(SpecAxis);
Range = TimeRange;
if (TimeRange(1) < 0)
    Range(1) = 0;
end
PlotSpectrogramInAxis(DirectoryName, FileName, FileType, SpecAxis, Range);
SpecAxisLimits = [TimeRange(1) TimeRange(2) 300 8000];
axis(SpecAxisLimits);
set(gca, 'YTick', []);
set(gca, 'XTickLabel', []);
ylabel('');

hold on;
plot([(TimeRange(1) + PreTime) (TimeRange(1) + PreTime)], [300 8000], 'k--');
plot([(TimeRange(2) - PostTime) (TimeRange(2) - PostTime)], [300 8000], 'k--');

plot([FileDur(1) FileDur(1)], [300 8000], 'b--');
plot([FileDur(2) FileDur(2)], [300 8000], 'b--');

cla(SpikeTrainAxis);
axes(SpikeTrainAxis);
plot(SpikeTimes, ones(size(SpikeTimes))*1, 'k+');
axis([TimeRange(1) TimeRange(2) 0.9 1.1]);
hold on;
plot([(TimeRange(1) + PreTime) (TimeRange(1) + PreTime)], [0.9 1.1], 'k--');
plot([(TimeRange(2) - PostTime) (TimeRange(2) - PostTime)], [0.9 1.1], 'k--');

plot([FileDur(1) FileDur(1)], [0.9 1.1], 'b--');
plot([FileDur(2) FileDur(2)], [0.9 1.1], 'b--');

set(gca, 'FontSize', 10);
set(gca, 'YTick', []);
xlabel('Time (sec)', 'FontSize', 12);


