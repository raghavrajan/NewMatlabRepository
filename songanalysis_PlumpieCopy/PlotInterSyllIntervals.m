function [] = PlotInterSyllIntervals(INR)

Index = 1;
for j = 1:length(INR.BoutDetails),
    for k = 1:length(INR.BoutDetails(j).labels) - 1,
        Labels{Index} = INR.BoutDetails(j).labels(k:k+1);
        Intervals(Index) = INR.BoutDetails(j).onsets(k+1) - INR.BoutDetails(j).offsets(k);
        Index = Index + 1;
    end
end

UniqueSylls = unique(cell2mat(Labels));

Index = 1;
IntervalIndices = ones(size(Intervals));
for i = 1:length(UniqueSylls),
    for j = 1:length(UniqueSylls),
        Indices = find(cellfun(@length, strfind(Labels, [UniqueSylls(i) UniqueSylls(j)])));
        if (length(Indices) >= 10)
            XLabelString{Index} = [UniqueSylls(i) UniqueSylls(j)];
            IntervalIndices(Indices) = Index;
            Index = Index + 1;
        end
    end
end

figure
boxplot(Intervals*1000, IntervalIndices);
hold on;
axis tight;
Temp = axis;
Temp(1:2) = [0.5 (length(XLabelString) + 0.5)];
Temp(4) = 500;
axis(Temp);

set(gca, 'XTick', [1:1:length(XLabelString)], 'XTickLabel', XLabelString);

% Now plot a line corresponding to the i to i duration
Index = find(cellfun(@length, strfind(XLabelString, 'ii')));
plot(Temp(1:2), median(Intervals(find(IntervalIndices == Index)))*ones(1,2)*1000, 'k--', 'LineWidth', 2);

