function [] = Harini_ProportionDirectedSongsAtEachLocation(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% ======= Proportion of directed songs at each location ===================
% Plots the proportion of songs at each location that are scored as
% directed. This is averaged across all experimental days for that
% condition.
% =========================================================================

Colours = 'rgbcmk';
Symbols = '+o<sd';

Conditions = [{'L0'} {'L1'} {'L2'} {'L3'} {'L4'} {'UN'}];

for i = 1:length(IndividualBirds),
    for j = 1:length(Conditions),
        eval([eval(['Conditions{', num2str(j), '}']), '(', num2str(i), ') = ', num2str(NaN), ';']);
    end
end
for i = 1:length(IndividualBirds),
    for j = 1:length(IndividualBirds(i).Conditions),
        % First find all the bouts that correspond to the given condition 
        ConditionIndices = find(IndividualBirds(i).BoutStatistics(:,end) == j);
        
        % Now within these, I have to see how many have video scoring data,
        % i.e. they don't have NA. - this is the total number of songs
        NumSongs(i,j) = length(ConditionIndices) - length(find(strcmp('NA', IndividualBirds(i).BoutCategorisation(ConditionIndices))));
        
        % Now find out how many are directed songs
        NumDirSongs(i,j) = length(find(strcmp('D', IndividualBirds(i).BoutCategorisation(ConditionIndices))));
        
        PropDirSongs = NumDirSongs(i,j) * 100 / NumSongs(i,j);
        
        eval([IndividualBirds(i).Conditions{j}, '(', num2str(i), ') = ', num2str(PropDirSongs), ';']);
    end
end

AllProportions = [L0(:) L1(:) L2(:) L3(:) L4(:) UN(:)];

for i = 1:length(IndividualBirds),
    NumConditions(i) = length(IndividualBirds(i).ConditionDistances);
end

[MaxNumConditions, MaxNumConditionsIndex] = max(NumConditions);

ConditionDistances = IndividualBirds(MaxNumConditionsIndex).ConditionDistances;

% Now plot the data
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 220 900 700]);
hold on;
for i = 1:size(AllProportions,2),
    AverageBar(i) = bar(ConditionDistances(i), mean(AllProportions(:,i), 1, 'default', 'omitnan'));
    set(AverageBar(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 2, 'BarWidth', 10);
    errorbar(ConditionDistances(i), mean(AllProportions(:,i), 1, 'default', 'omitnan'), std(AllProportions(:,i), 0, 1, 'omitnan')/sqrt(length(find(~isnan(AllProportions(:,i))))), 'ko-', 'LineWidth', 2, 'MarkerSize', 10);
end

for i = 1:size(AllProportions, 1),
    plot(ConditionDistances, AllProportions(i,:), [Colours(mod(i, length(Colours)) + 1), Symbols(mod(i, length(Symbols)) + 1), '-'], 'LineWidth', 0.5, 'MarkerSize', 8);
end

axis tight;
Temp = axis;
Temp = [(ConditionDistances(1)-15) (ConditionDistances(end)+15) 0 min(1.05*Temp(4), 100)];
axis(Temp);
set(gca, 'XTick', ConditionDistances, 'XTickLabel', [mat2cell(ConditionDistances(1:5), 1, ones(5,1)), {'UN'}]);
xlabel('Distance of female');
ylabel('Proportion of songs that were scored as DIRECTED');



