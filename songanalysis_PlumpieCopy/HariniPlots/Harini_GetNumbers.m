function [DirBout_Stats, UndirBout_Stats, DirUndirBout_Stats] = Harini_GetNumbers(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% Script for getting the number of bouts, etc. across all birds.

for i = 1:length(IndividualBirds),
    disp(['Bird #', num2str(i)]);
    ConditionColumnIndex = find(cellfun(@length,strfind(IndividualBirds(i).BoutStatisticsColumnNames, 'Condition')));
    BoutIndexColumnIndex = find(cellfun(@length,strfind(IndividualBirds(i).BoutStatisticsColumnNames, 'BoutIndex')));
    
    INLabels = IndividualBirds(i).SortedBirdParameters(1).INLabels;
    MotifLabels = IndividualBirds(i).SortedBirdParameters(1).MotifLabels;
    CommonMotifs = IndividualBirds(i).SortedBirdParameters(1).CommonMotifs;
    
    ValidSongBouts = find((IndividualBirds(i).Bouts(:,7) == 1) & (IndividualBirds(i).Bouts(:,8) > 0) & (IndividualBirds(i).Bouts(:,9) > 1));
    disp(['% of song bouts that have enough data in front and at the back = ', num2str(length(ValidSongBouts)), '/', num2str(length(find(IndividualBirds(i).Bouts(:,7) == 1))), '(', num2str(100 * length(ValidSongBouts)/length(find(IndividualBirds(i).Bouts(:,7) == 1))), '%)']);
    for j = min(IndividualBirds(i).BoutStatistics(:,ConditionColumnIndex)):max(IndividualBirds(i).BoutStatistics(:,ConditionColumnIndex)),
        Bouts{i}{j} = IndividualBirds(i).BoutStatistics(find(IndividualBirds(i).BoutStatistics(:,ConditionColumnIndex) == j), BoutIndexColumnIndex);
        Bouts{i}{j} = intersect(Bouts{i}{j}, ValidSongBouts);
        
        DirBouts{i}{j} = intersect(IndividualBirds(i).BoutStatistics(strmatch('D', IndividualBirds(i).BoutCategorisation, 'exact'), BoutIndexColumnIndex), Bouts{i}{j});
        DirUndirBouts{i}{j} = intersect(IndividualBirds(i).BoutStatistics(strmatch('DUN', IndividualBirds(i).BoutCategorisation, 'exact'), BoutIndexColumnIndex), Bouts{i}{j});
        if (j ~= 6)
            UndirBouts{i}{j} = intersect(IndividualBirds(i).BoutStatistics(strmatch('UN', IndividualBirds(i).BoutCategorisation, 'exact'), BoutIndexColumnIndex), Bouts{i}{j});
        else
            UndirBouts{i}{j} = Bouts{i}{j}; % since these bouts are anyway only in the absence of the female.
        end
        TotalBoutNos(i,j) = length(DirBouts{i}{j}) + length(DirUndirBouts{i}{j}) + length(UndirBouts{i}{j});
    end
    DirBoutNos(i,:) = cellfun(@length, DirBouts{i});
    UnDirBoutNos(i,:) = cellfun(@length, UndirBouts{i});
end

%TotalBoutNos(find(TotalBoutNos == 0)) = NaN;
NonNaNBouts = find(~isnan(TotalBoutNos));
disp(['Median # of bouts / distance / bird = ', num2str(median(TotalBoutNos(NonNaNBouts))), '; range = ', num2str(min(TotalBoutNos(NonNaNBouts))), ' - ', num2str(max(TotalBoutNos(NonNaNBouts)))]);

MinTrialNo = 3;
% Now display number of distances for each bird with more than 3 
DirNosMoreThanMin = DirBoutNos >= MinTrialNo;
UnDirNosMoreThanMin = UnDirBoutNos >= MinTrialNo;

% Exclude Bird #6 and Bird #14 because they don't have L0 or have only L0
% and UN.

DirDistanceMoreThanMin = sum(DirNosMoreThanMin, 2);
UnDirDistanceMoreThanMin = sum(UnDirNosMoreThanMin(:,1:5), 2);

DirDistanceMoreThanMin([6, 14]) = [];
UnDirDistanceMoreThanMin([6, 14]) = [];
disp(['Dir songs: >= ', num2str(MinTrialNo), ' song bouts: Median # of distances / bird = ', num2str(median(DirDistanceMoreThanMin)), ', range = ', num2str(min(DirDistanceMoreThanMin)), ' - ', num2str(max(DirDistanceMoreThanMin))]);
disp(['UnDir songs: >= ', num2str(MinTrialNo), ' song bouts: Median # of distances / bird = ', num2str(median(UnDirDistanceMoreThanMin)), ', range = ', num2str(min(UnDirDistanceMoreThanMin)), ' - ', num2str(max(UnDirDistanceMoreThanMin))]);