function [] = CSOverlayAcousticFeaturesPlot(CSData)

FeatureAxisLabels = [{'Duration'} {'Duration (sec)'}; {'LogAmplitude'} {'Amplitude (dB)'}; {'Entropy'} {'Entropy'}; {'PitchGoodness'} {'Pitch Goodness'}; {'FrequencyModulation'} {'Frequency Modulation'}; {'AmplitudeModulation'} {'Amplitude Modulation (dB/sec)'}; {'FundamentalFrequency'} {'Frequency (Hz)'}; {'MeanFrequency'} {'Mean frequency (Hz)'}];

XFeature = inputdlg('Choose the feature that you want to plot on the x-axis', 'X-axis Feature');
YFeature = inputdlg('Choose the feature that you want to plot on the y-axis', 'Y-axis Feature');

DaysToPlot = inputdlg('Choose the days that you want plotted: Just type the numbers separated by commas', 'Days to be plotted');

CommaIndices = find(DaysToPlot{1} == ',');

if (~isempty(CommaIndices))
    ActualDaysToPlot(1) = str2double(DaysToPlot{1}(1:CommaIndices(1)-1));
    for i = 1:length(CommaIndices),
        if (i == length(CommaIndices))
            ActualDaysToPlot(end+1) = str2double(DaysToPlot{1}(CommaIndices(i)+1:end));
        else
            ActualDaysToPlot(end+1) = str2double(DaysToPlot{1}(CommaIndices(i)+1:CommaIndices(i+1)-1));
        end
    end
else
    ActualDaysToPlot(1) = str2double(DaysToPlot{1});
end

OverlayFigure = figure;
Legend = [];
for i = ActualDaysToPlot,
    XFeatureIndex = strmatch(XFeature{1}, CSData.Data{i}.ToBeUsedFeatures, 'exact');
    YFeatureIndex = strmatch(YFeature{1}, CSData.Data{i}.ToBeUsedFeatures, 'exact');
    Color = uisetcolor([], ['Day #', num2str(i), ' Colour']);
    Temp = inputdlg(['Enter the legend for Day #', num2str(i)], 'Legend');
    Legend{end+1} = Temp{1};
    figure(OverlayFigure);
    plot(CSData.AllFeats{i}(:, XFeatureIndex), CSData.AllFeats{i}(:, YFeatureIndex), 'ko', 'Color', Color);
    hold on;
end

figure(OverlayFigure);
set(gcf, 'Color', 'w');
set(gca, 'FontSize', 16);
XAxisLabelIndex = strmatch(XFeature{1}, FeatureAxisLabels(:,1), 'exact');
xlabel(FeatureAxisLabels{XAxisLabelIndex, 2}, 'FontSize', 16);

YAxisLabelIndex = strmatch(YFeature{1}, FeatureAxisLabels(:,1), 'exact');
ylabel(FeatureAxisLabels{YAxisLabelIndex, 2}, 'FontSize', 16);

legend(Legend);

disp('Finished plotting data');