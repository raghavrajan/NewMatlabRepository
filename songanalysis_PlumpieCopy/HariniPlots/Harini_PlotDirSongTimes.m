function [] = Harini_PlotDirSongTimes(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% ======= Proportion of directed songs at each location ===================
% Plots the proportion of songs at each location that are scored as
% directed. This is averaged across all experimental days for that
% condition.
% =========================================================================

Colours = 'rgbcmk';
Symbols = '+o<sd';

for i = 1:length(IndividualBirds),
    figure;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [680 220 900 780]);
    hold on;
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
        
        plot(ones(size(DirSongIndices))*j, IndividualBirds(i).BoutStatisticsTimeRelativeToStart(DirSongIndices)*60, 'k+');
    end
    
    axis tight;
    Temp = axis;
    Temp = [0.5 j+0.5 0 60];
    axis(Temp);
    set(gca, 'XTick', 1:1:j, 'XTickLabel', IndividualBirds(i).Conditions);
    xlabel('Condition');
    ylabel('Time of occurence of DIRECTED song relative to start of the session (mins)');
    title(BirdNames{i});
end

