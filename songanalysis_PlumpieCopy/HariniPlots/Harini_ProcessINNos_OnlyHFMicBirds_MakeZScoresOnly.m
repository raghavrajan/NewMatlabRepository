function [DirBout_Stats, UndirBout_Stats, DirUndirBout_Stats] = Harini_ProcessINNos_OnlyHFMicBirds_MakeZScoresOnly(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% ====== IN Number processing =============================================
% THis is to check all aspects related to IN number in all the birds.

% =========================================================================

FinalFigureDir = '/home/raghav/StudentRelated/Harini/Manuscript/results/OutliersRemoved';
MinTrialNo = 3;

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
            [DirBout_Stats{i}{j}] = Harini_CalculateBoutINMotifNumbers_Duration(IndividualBirds(i), DirBouts{i}{j}, MotifLabels, INLabels, CommonMotifs);
        else
            DirBout_Stats{i}{j} = [];
        end
        if (~isempty(DirUndirBouts{i}{j}))
            [DirUndirBout_Stats{i}{j}] = Harini_CalculateBoutINMotifNumbers_Duration(IndividualBirds(i), DirUndirBouts{i}{j}, MotifLabels, INLabels, CommonMotifs);
        else
            DirUndirBout_Stats{i}{j} = [];
        end
        if (~isempty(UndirBouts{i}{j}))
            [UndirBout_Stats{i}{j}] = Harini_CalculateBoutINMotifNumbers_Duration(IndividualBirds(i), UndirBouts{i}{j}, MotifLabels, INLabels, CommonMotifs);
        else
            UndirBout_Stats{i}{j} = [];
        end
    end
end


Distances = [0 20 60 110 165 200];

MinTrialNo = 10;
StructString = 'LogAmplitude';
LabelString = 'Log Amplitude (dB)';

Harini_AmplitudeZScores_Plots_OnlyDir(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames);

disp('Finished calculating bout stats');