function [] = Harini_PlotINNumber_DirOnlySongs(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% ======= IN number for only the songs classified as DIRECTED =============
% Plots the number of INs for only songs classified as directed at each of
% the locations - averaged over multiple days of experiments
% =========================================================================

Colours = 'rgbcmk';
Symbols = '+o<sd';


for i = 1:length(IndividualBirds),
    ColumnIndex = find(strcmp('BoutLength', IndividualBirds(i).BoutStatisticsColumnNames));
    % Now plot the distributions
    figure;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [680 220 700 500]);
    hold on;
    PlotLegend{i} = [];
    for j = 1:length(IndividualBirds(i).Conditions),
        % First find all the bouts that correspond to the given condition 
        ConditionIndices = find(IndividualBirds(i).BoutStatistics(:,end) == j);
        
        % Now within this set, find the ones that are actually directed
        % This is only for all the directed conditions - for undirected
        % condition, just plot everything
        if (isempty(strfind(IndividualBirds(i).Conditions{j}, 'UN')))
            DirSongIndices = ConditionIndices(find(strcmp('D', IndividualBirds(i).BoutCategorisation(ConditionIndices))));
        else
            DirSongIndices = ConditionIndices;
        end
        if (length(DirSongIndices) >= 3)
            plot(0:500:max(IndividualBirds(i).BoutStatistics(:, ColumnIndex)), histc(IndividualBirds(i).BoutStatistics(DirSongIndices, ColumnIndex), 0:500:max(IndividualBirds(i).BoutStatistics(:,ColumnIndex)))/length(DirSongIndices), [Colours(mod(j, length(Colours)) + 1), 'o-'], 'LineWidth', 2);
            PlotLegend{i}{end+1} = IndividualBirds(i).Conditions{j};
        end
    end
    
    for j = 1:length(IndividualBirds(i).Conditions),
        ConditionIndices = find(IndividualBirds(i).BoutStatistics(:,end) == j);
        if (isempty(strfind(IndividualBirds(i).Conditions{j}, 'UN')))
            DirSongIndices = ConditionIndices(find(strcmp('D', IndividualBirds(i).BoutCategorisation(ConditionIndices))));
        else
            DirSongIndices = ConditionIndices;
        end
        if (length(DirSongIndices) >= 3)
            plot(mean(IndividualBirds(i).BoutStatistics(DirSongIndices, ColumnIndex)), 1 + j*0.05, [Colours(mod(j, length(Colours)) + 1), 's'], 'MarkerSize', 8, 'LineWidth', 1.5);
            plot([(mean(IndividualBirds(i).BoutStatistics(DirSongIndices, ColumnIndex)) - std(IndividualBirds(i).BoutStatistics(DirSongIndices, ColumnIndex))) (mean(IndividualBirds(i).BoutStatistics(DirSongIndices, ColumnIndex)) + std(IndividualBirds(i).BoutStatistics(DirSongIndices, ColumnIndex)))], [(1 + j*0.05) (1 + j*0.05)], Colours(mod(j, length(Colours)) + 1), 'LineWidth', 1.5);
        end
    end
    
    axis tight;
    Temp = axis;
    Temp = [0 (max(IndividualBirds(i).BoutStatistics(:,ColumnIndex)) + 0.5) 0 1.35];
    axis(Temp);
    legend(PlotLegend{i});
    xlabel('Bout length (msec)');
    ylabel('Fraction of bouts');
    title(BirdNames{i});
end


