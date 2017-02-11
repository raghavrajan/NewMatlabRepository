function [] = CSAnalyseSyllableAcousticFeatures(CSData)

UniqueLabels = [];
for i = 1:CSData.NoofDays,
    UniqueLabels = union(UniqueLabels, unique(CSData.AllLabels{i}));
end
UniqueLabels = UniqueLabels(:)';
UniqueLabels = UniqueLabels(regexp(UniqueLabels, '[a-p 0 r-z]'));

for i = 1:length(CSData.Data{1}.ToBeUsedFeatures),
    figure(i);
    hold on;
    MeanFeatData = [];
    STDFeatData = [];
    SEMFeatData = [];
    for j = 1:length(UniqueLabels),
        for k = 1:CSData.NoofDays,
            Indices = find(CSData.AllLabels{k} == UniqueLabels(j));
            if (~isempty(Indices))
                MeanFeatData(j, k) = mean(CSData.AllFeats{k}(Indices, i));
                STDFeatData(j, k) = std(CSData.AllFeats{k}(Indices, i));
                SEMFeatData(j, k) = STDFeatData(j, k)/sqrt(length(Indices));
            else
                MeanFeatData(j, k) = NaN;
                STDFeatData(j, k) = NaN;
                SEMFeatData(j, k) = NaN;
            end
        end
    end
    BarPlotHandle = bar(MeanFeatData);
    hold on;
    colormap('gray');
    for j = 1:size(MeanFeatData, 2),
        XVal = get(get(BarPlotHandle(j), 'children'), 'xData');
        errorbar(mean(XVal,1), MeanFeatData(:,j), STDFeatData(:,j), 'k.', 'MarkerSize', 2);
    end
    axis tight;
    set(gca, 'XTick', 1:1:length(UniqueLabels), 'XTickLabel', mat2cell(UniqueLabels, 1, ones(size(UniqueLabels))), 'FontSize', 16);
    title(CSData.Data{1}.ToBeUsedFeatures{i}, 'FontSize', 16);
end
        
disp('Finished plotting syllable acoustic features');