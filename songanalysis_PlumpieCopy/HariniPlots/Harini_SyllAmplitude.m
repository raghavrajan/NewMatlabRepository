function [DirBout_Stats, UndirBout_Stats, DirUndirBout_Stats] = Harini_SyllAmplitude(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% ====== Final Analysis Script ============================================
% This is to analyze and plot the following parameters
% 1. # of INs
% 2. # of complete motifs / bout
% 3. All motif duration
% 4. CV of FF for each individual syllable
% 5. Log amplitude for each individual syllable
% =========================================================================

FinalFigureDir = '/home/raghav/StudentRelated/Harini/Manuscript/results';

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
        % Now for each of these bouts, get the # of INs, the last set of INs
        % with < 500ms gap and the last set of INs independent of the gap
        % In addition, for each bout, also get the total number of complete
        % motifs, partial motifs, all motif durations, first motif duration
        % Also, get the date for that particular dataset
    
        if (~isempty(DirBouts{i}{j}))
            [DirBout_Stats{i}{j}] = Harini_CalculateBoutINMotifNumbers_Duration_AmpKao(IndividualBirds(i), DirBouts{i}{j}, MotifLabels, INLabels, CommonMotifs);
        else
            DirBout_Stats{i}{j} = [];
        end
        if (~isempty(DirUndirBouts{i}{j}))
            [DirUndirBout_Stats{i}{j}] = Harini_CalculateBoutINMotifNumbers_Duration_AmpKao(IndividualBirds(i), DirUndirBouts{i}{j}, MotifLabels, INLabels, CommonMotifs);
        else
            DirUndirBout_Stats{i}{j} = [];
        end
        if (~isempty(UndirBouts{i}{j}))
            [UndirBout_Stats{i}{j}] = Harini_CalculateBoutINMotifNumbers_Duration_AmpKao(IndividualBirds(i), UndirBouts{i}{j}, MotifLabels, INLabels, CommonMotifs);
        else
            UndirBout_Stats{i}{j} = [];
        end
    end
end

% First get number of bouts for all types
for i = 1:length(DirBout_Stats),
    DirBoutNos(i,:) = cellfun(@length, DirBouts{i});
    DirUndirBoutNos(i,:) = cellfun(@length, DirUndirBouts{i});
    UndirBoutNos(i,:) = cellfun(@length, UndirBouts{i});
end

Distances = [0 20 60 110 165 195];
DistanceLabelString = {'0' '20' '60' '110' '165' 'UN'};

% Now for each bird plot the # of INs, # of motifs/bout and motif duration
% with the corresponding fits for individual birds

MinTrialNo = 3;

% Now find all birds that have atleast with # of trials for L0 > MinTrialNo
% and atleast 3 distances (including) L0 wth more than the MinTrialNo

for i = 1:size(DirBoutNos,1),
    NumDistances(i) = length(find(DirBoutNos(i,:) >= MinTrialNo));
end

% Use only birds that have data for at least 3 distances
BirdsToUse = find((NumDistances(:) >= 3));
disp(['Number of birds with >= ', num2str(MinTrialNo), ' directed song bouts at 3 distances = ', num2str(length(BirdsToUse)), ' / ', num2str(size(DirBoutNos, 1)), ' birds']);
disp(['Median # of distances / bird = ', num2str(median(NumDistances(BirdsToUse))), '; range = ', num2str(min(NumDistances(BirdsToUse))), ' - ', num2str(max(NumDistances(BirdsToUse)))]);

NumBoutsPerDistance = [];
for i = BirdsToUse(:)',
    ValidDistances = find(DirBoutNos(i,:) >= MinTrialNo);
    NumBoutsPerDistance = [NumBoutsPerDistance(:); DirBoutNos(i, ValidDistances)'];
end
disp(['Median # of directed song bouts / distance = ', num2str(median(NumBoutsPerDistance)), '; range = ', num2str(min(NumBoutsPerDistance)), ' - ', num2str(max(NumBoutsPerDistance))]);


Harini_FeatureFits_Amplitude(DirBout_Stats(BirdsToUse), DirUndirBout_Stats(BirdsToUse), UndirBout_Stats(BirdsToUse), Distances, BirdNames(BirdsToUse));

disp('Finished calculating bout stats');