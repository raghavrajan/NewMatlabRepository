function [] = LSINA_PlotProportionDiversitySyllables(DataStruct, PlotString, FileName)

% This function plots the diversity and proportion of syllables given the
% input data
% Written primarily for plotting the proportion and diversity of syllables
% found at the beginning of a bout, beginning of a motif and the syllable
% just before the motif.

BoutStartSylls = zeros(length(DataStruct), max([DataStruct.NumDiffSylls]));
for i = 1:length(DataStruct),
    BoutStartSylls(i,1:DataStruct(i).NumDiffSylls) = DataStruct(i).ProportionSylls;
end
figure;
hold on;
bar(BoutStartSylls, 'stacked');
for i = 1:length(DataStruct),
    XAxis_LabelStrings{i} = [DataStruct(i).BirdName, '.', DataStruct(i).DataLabel];
    for j = 1:length(DataStruct(i).DiffSylls),
        if (j == 1)
            text(i, DataStruct(i).ProportionSylls(j)/2, DataStruct(i).DiffSylls(j), 'FontSize', 12, 'FontWeight', 'bold');
        else
            text(i, sum(DataStruct(i).ProportionSylls(1:j-1)) + DataStruct(i).ProportionSylls(j)/2, DataStruct(i).DiffSylls(j), 'FontSize', 12, 'FontWeight', 'bold');
        end
    end
end
PlotAxis = [0 (length(DataStruct) + 1) 0 1.1];

LSINA_FixPlotLabels(gcf, 10, 'Fig_Position', [94 195 1100 350], 'Fig_Title', ['Diversity and proportion of ', PlotString, ' syllables'], 'YAxis_Label', 'Fraction', 'Fig_XTickLabels', XAxis_LabelStrings, 'Fig_XTickLabel_Rotation', 45, 'Fig_Axis', PlotAxis, 'Fig_Save', FileName);