function [FeatureVect, GroupVect] = INAPlotFeatures(Feature1, Feature2, FeatureName)

FeatureVect = [];
GroupVect = [];
for i = 1:length(Feature1),
    FeatLen(i) = length(Feature1{i});
end

MaxFeatLen = find(FeatLen >= 5, 1, 'last');

for i = 1:MaxFeatLen,
    Bar(i) = bar(i, median(Feature1{i}));
    set(Bar(i), 'FaceColor', 'w', 'EdgeColor', 'k');
    hold on;
    plot(ones(size(Feature1{i}))*i, Feature1{i}, 'k+', 'MarkerSize', 2);
    FeatureVect = [FeatureVect; Feature1{i}];
    GroupVect = [GroupVect; ones(size(Feature1{i}))*i];
end
if (~isempty(Feature2))
    Bar(6) = bar(6, mean(Feature2));
    set(Bar(6), 'FaceColor', 'w', 'EdgeColor', 'b');
    plot(ones(size(Feature2))*6, Feature2, 'k+', 'MarkerSize', 2);
    FeatureVect = [FeatureVect; Feature2];
    GroupVect = [GroupVect; ones(size(Feature2))*6];
end
axis tight
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
ylabel(FeatureName, 'FontSize', 14, 'FontWeight', 'bold');

[p, anovatab, stats] = anova1(FeatureVect, GroupVect);
title(FeatureName, 'FontSize', 14, 'FontWeight', 'bold');
multcompare(stats);
title(FeatureName, 'FontSize', 14, 'FontWeight', 'bold');

