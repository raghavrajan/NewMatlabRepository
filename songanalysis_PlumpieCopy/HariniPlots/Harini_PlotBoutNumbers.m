function [] = Harini_PlotBoutNumbers(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% ======= Display the total number of bouts for all conditions ============
% 1. All bouts
% 2. Directed only bouts
% 3. Undirected bouts
% =========================================================================

Colours = 'rgbcmk';
Symbols = '+o<sd';

Fid = fopen('AllBirds_BoutNumbers.txt', 'w');
fprintf(Fid, 'Bird #\tBird Name\tCondition\tTotal # of bouts\t# of bouts with video scoring\t# of Directed bouts\t# of Undirected bouts\n');
for i = 1:length(IndividualBirds),
    ConditionColumnIndex = find(strcmp('Condition', IndividualBirds(i).BoutStatisticsColumnNames));
    for j = 1:length(IndividualBirds(i).Conditions),
        % First find all the bouts that correspond to the given condition 
        ConditionIndices = find(IndividualBirds(i).BoutStatistics(:,ConditionColumnIndex) == j);
        % Find all bouts with video scoring data - anything but 'NA'
        if (isempty(strfind(IndividualBirds(i).Conditions{j}, 'UN')))
            VideoScoredBoutIndices = ConditionIndices(find(~strcmp('NA', IndividualBirds(i).BoutCategorisation(ConditionIndices))));
        else
            VideoScoredBoutIndices = ConditionIndices;
        end
        % Find all bouts with 'D' label - i.e. directed songs as per video
        % scoring
        if (isempty(strfind(IndividualBirds(i).Conditions{j}, 'UN')))
            DirSongIndices = ConditionIndices(find(strcmp('D', IndividualBirds(i).BoutCategorisation(ConditionIndices))));
        else
            DirSongIndices = ConditionIndices;
        end
        % Find all bouts with 'UN' label - i.e. undirected songs as per
        % video scoring
        if (isempty(strfind(IndividualBirds(i).Conditions{j}, 'UN')))
            UnDirSongIndices = ConditionIndices(find(strcmp('UN', IndividualBirds(i).BoutCategorisation(ConditionIndices))));
        else
            UnDirSongIndices = ConditionIndices;
        end
        fprintf(Fid, '%d\t%s\t%s\t%d\t%d\t%d\t%d\n', i, BirdNames{i}, IndividualBirds(i).Conditions{j}, length(ConditionIndices), length(VideoScoredBoutIndices), length(DirSongIndices), length(UnDirSongIndices));
    end
end
fclose(Fid);
disp('Finished writing number of bout info to file');
    