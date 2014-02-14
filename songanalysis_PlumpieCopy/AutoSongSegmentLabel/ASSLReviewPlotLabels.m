function [] = ASSLReviewPlotLabels(Onsets, Labels, LabelAxis, LabelAxisLimits, CurrentLabelIndex)

Onsets = Onsets/1000;

cla(LabelAxis);
axes(LabelAxis);
for i = 1:length(Labels),
    if (i ~= CurrentLabelIndex)
        text(Onsets(i), 1, Labels(i), 'Color', 'b');
    else
        text(Onsets(i), 1, Labels(i), 'Color', 'r');
    end
end
axis(LabelAxisLimits);
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
set(gca, 'XTick', []);
set(gca, 'YTick', []);