function [] = CSAnalyseINLikeSyllables(CSData)

SyllablesToAnalyse = inputdlg('Input the syllables (starting from the beginning of the bout with 1) that you want to analyse (use comma to separate multiple syllable nos.)', 'IN-like syllables');

CommaIndices = find(SyllablesToAnalyse{1} == ',');

if (~isempty(CommaIndices))
    ActualSyllIndices(1) = str2double(SyllablesToAnalyse{1}(1:CommaIndices(1)-1));
    for i = 1:length(CommaIndices),
        if (i == length(CommaIndices))
            ActualSyllIndices(end+1) = str2double(SyllablesToAnalyse{1}(CommaIndices(i)+1:end));
        else
            ActualSyllIndices(end+1) = str2double(SyllablesToAnalyse{1}(CommaIndices(i)+1:CommaIndices(i+1)-1));
        end
    end
else
    ActualSyllIndices(1) = str2double(SyllablesToAnalyse{1});
end

for i = 1:CSData.NoofDays,
    Feats{i} = [];
    Intervals{i} = [];
    BoutStarts = find(CSData.AllLabels{i} == 'Q');
    BoutEnds = find(CSData.AllLabels{i} == 'q');
    for j = 1:length(BoutStarts),
       MotifSylls = regexp(CSData.AllLabels{i}(BoutStarts(j):BoutEnds(j)), ['[', CSData.MotifSyllLabels, ']']);
       if (isempty(MotifSylls))
           continue;
       end
       Feats{i}(end+1:(end+length(ActualSyllIndices)), :) = CSData.AllFeats{i}(BoutStarts(j) + ActualSyllIndices, :);
       if (length(ActualSyllIndices) > 1)
           Intervals{i}(end+1, :) = (CSData.AllOnsets{i}(BoutStarts(j) + ActualSyllIndices(2:end)) - CSData.AllOffsets{i}(BoutStarts(j) + ActualSyllIndices(1:end-1)))';
       end
    end
end


FeatureAxisLabels = [{'Duration'} {'Duration (sec)'}; {'LogAmplitude'} {'Amplitude (dB)'}; {'Entropy'} {'Entropy'}; {'PitchGoodness'} {'Pitch Goodness'}; {'FrequencyModulation'} {'Frequency Modulation'}; {'AmplitudeModulation'} {'Amplitude Modulation (dB/sec)'}; {'FundamentalFrequency'} {'Frequency (Hz)'}; {'MeanFrequency'} {'Mean frequency (Hz)'}];

XFeature = inputdlg('Choose the feature that you want to plot on the x-axis', 'X-axis Feature');
YFeature = inputdlg('Choose the feature that you want to plot on the y-axis', 'Y-axis Feature');

FeatPlotFigure = figure;
Legend = [];
for i = 1:CSData.NoofDays,
    XFeatureIndex = strmatch(XFeature{1}, CSData.Data{i}.ToBeUsedFeatures, 'exact');
    YFeatureIndex = strmatch(YFeature{1}, CSData.Data{i}.ToBeUsedFeatures, 'exact');
    Color{i} = uisetcolor([], ['Day #', num2str(i), ' Colour']);
    Temp = inputdlg(['Enter the legend for Day #', num2str(i)], 'Legend');
    Legend{end+1} = Temp{1};
    figure(FeatPlotFigure);
    plot(Feats{i}(:, XFeatureIndex), Feats{i}(:, YFeatureIndex), 'ko', 'Color', Color{i});
    hold on;
end

figure(FeatPlotFigure);
set(gcf, 'Color', 'w');
set(gca, 'FontSize', 16);
XAxisLabelIndex = strmatch(XFeature{1}, FeatureAxisLabels(:,1), 'exact');
xlabel(FeatureAxisLabels{XAxisLabelIndex, 2}, 'FontSize', 16);

YAxisLabelIndex = strmatch(YFeature{1}, FeatureAxisLabels(:,1), 'exact');
ylabel(FeatureAxisLabels{YAxisLabelIndex, 2}, 'FontSize', 16);

legend(Legend);

IntervalFigure = figure;
for i = 1:CSData.NoofDays,
    Intervals{i} = sort(Intervals{i});
    %plot(median(Intervals{i}), 'ko', 'Color', Color{i});
    hold on;
    errorbar(1:1:size(Intervals{i},2), median(Intervals{i}), median(Intervals{i}) - Intervals{i}(round(0.25*size(Intervals{i},1)),:), Intervals{i}(round(0.75*size(Intervals{i},1)),:) - median(Intervals{i}), 'ko-', 'Color', Color{i});
end
set(gcf, 'Color', 'w');
set(gca, 'FontSize', 16);
for i = 1:size(Intervals{i},2),
    XAxisLabel{i} = [num2str(ActualSyllIndices(i)), '-', num2str(ActualSyllIndices(i+1))];
end
set(gca, 'XTick', 1:1:size(Intervals{1},2), 'XTickLabel', XAxisLabel);
ylabel('Median interval (msec)', 'FontSize', 16);
legend(Legend);

disp('Finished Analysis');