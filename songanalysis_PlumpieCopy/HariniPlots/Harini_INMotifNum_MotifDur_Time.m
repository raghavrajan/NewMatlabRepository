function [DirBout_Stats, UndirBout_Stats, DirUndirBout_Stats] = Harini_INMotifNum_MotifDur_Time(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% ====== Final Analysis Script ============================================
% This is to analyze and plot the following parameters
% 1. # of INs
% 2. # of complete motifs / bout
% 3. All motif duration
% 4. CV of FF for each individual syllable
% 5. Log amplitude for each individual syllable
% =========================================================================

DirSongEndTime = 10; % in minutes
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
        AllBoutTimes = [];
        AllBouts = [];
        
        Bouts{i}{j} = IndividualBirds(i).BoutStatistics(find(IndividualBirds(i).BoutStatistics(:,ConditionColumnIndex) == j), BoutIndexColumnIndex);
        Bouts{i}{j} = intersect(Bouts{i}{j}, ValidSongBouts);
        
        DirBouts{i}{j} = intersect(IndividualBirds(i).BoutStatistics(strmatch('D', IndividualBirds(i).BoutCategorisation, 'exact'), BoutIndexColumnIndex), Bouts{i}{j});
        DirUndirBouts{i}{j} = intersect(IndividualBirds(i).BoutStatistics(strmatch('DUN', IndividualBirds(i).BoutCategorisation, 'exact'), BoutIndexColumnIndex), Bouts{i}{j});
        
        if (j ~= 6)
            UndirBouts{i}{j} = intersect(IndividualBirds(i).BoutStatistics(strmatch('UN', IndividualBirds(i).BoutCategorisation, 'exact'), BoutIndexColumnIndex), Bouts{i}{j});
        else
            UndirBouts{i}{j} = Bouts{i}{j}; % since these bouts are anyway only in the absence of the female.
        end
        
        if (~isempty(DirBouts{i}{j}))
            DirBoutTimes{i}{j} = IndividualBirds(i).BoutTimesRelativeToStart(DirBouts{i}{j}) * 60; % in minutes
            AllBoutTimes = [AllBoutTimes; DirBoutTimes{i}{j}(:)];
            AllBouts = [AllBouts; DirBouts{i}{j}(:)];
        end
        if (~isempty(DirUndirBouts{i}{j}))
            DirUndirBoutTimes{i}{j} = IndividualBirds(i).BoutTimesRelativeToStart(DirUndirBouts{i}{j}) * 60; % in minutes
            AllBouts = [AllBouts; DirUndirBouts{i}{j}(:)];
            AllBoutTimes = [AllBoutTimes; DirUndirBoutTimes{i}{j}(:)];
        end
        if (~isempty(UndirBouts{i}{j}))
            UndirBoutTimes{i}{j} = IndividualBirds(i).BoutTimesRelativeToStart(UndirBouts{i}{j}) * 60; % in minutes
            AllBouts = [AllBouts; UndirBouts{i}{j}(:)];
            AllBoutTimes = [AllBoutTimes; UndirBoutTimes{i}{j}(:)];
        end
        
        % Based on bout times choose dir and undir bouts 
        DirBouts{i}{j} = AllBouts(find(AllBoutTimes <= DirSongEndTime));
 
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

% Now plot the percentage of dir songs at each distance

PercentageDirBoutNos = 100 * DirBoutNos./(DirBoutNos + UndirBoutNos);
DistanceMatrix = repmat(Distances(:)', size(PercentageDirBoutNos,1), 1);
NonNaNPercentages = find(~isnan(PercentageDirBoutNos));
[Corr_R, Corr_P] = corr(DistanceMatrix(NonNaNPercentages), PercentageDirBoutNos(NonNaNPercentages));

p = panel();
p.pack({1});
p(1).select();
hold on;
plot(repmat(Distances, size(PercentageDirBoutNos, 1), 1)', PercentageDirBoutNos', 'ko-', 'Color', [0.75 0.75 0.75]);
for i = 1:size(PercentageDirBoutNos, 2),
    Indices = find(~isnan(PercentageDirBoutNos(:,i)));
    errorbar(Distances(i), mean(PercentageDirBoutNos(Indices,i)), std(PercentageDirBoutNos(Indices,i))/sqrt(length(Indices)), 'ks', 'MarkerSize', 8, 'LineWidth', 2);
end
text(177, 0, '//', 'FontSize', 16)
axis([-5 Distances(end)+5 0 105]);
set(gca, 'Box', 'off')
set(gcf, 'Color', 'w')

xlabel({'Distance from the female (cm)'});
set(gca, 'XTick', Distances, 'XTickLabel', DistanceLabelString);
ylabel('Proportion of directed songs (%)');
text(Distances(round(length(Distances)/2)+1), 95, {['r = ', num2str(Corr_R)]; ['p = ', num2str(Corr_P)]});
p.fontsize = 12;
% p.fontname = 'Times';
% p.margintop = 10;
p.marginleft = 20;
% p.marginright = 20;

set(gcf, 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position', [680 621 425 350]);
print(fullfile(FinalFigureDir, ['PercentageofDirSongsWithDistance_DBasedOnTime.eps']), '-depsc2', '-r600');
print(fullfile(FinalFigureDir, ['PercentageofDirSongsWithDistance_DBasedOnTime.png']), '-dpng', '-r600');
close all;

% Now for each bird plot the # of INs, # of motifs/bout and motif duration
% with the corresponding fits for individual birds

MinTrialNo = 3;

% Now find all birds that have atleast with # of trials for L0 > MinTrialNo
% and atleast 3 distances (including) L0 wth more than the MinTrialNo

for i = 1:size(DirBoutNos,1),
    NumDistances(i) = length(find(DirBoutNos(i,:) >= MinTrialNo));
end

% First plot differences in feature values for L0 dir and Undir (alone) -
% use all birds that have >= 3 bouts for both categories

BirdsToUse = find((DirBoutNos(:,1) >= MinTrialNo) & (UndirBoutNos(:,end) >= MinTrialNo));
disp(['Number of birds with >= ', num2str(MinTrialNo), ' directed song bouts at L0 and >= ', num2str(MinTrialNo), ' undirected song bouts = ', num2str(length(BirdsToUse)), ' / ', num2str(size(DirBoutNos, 1)), ' birds']);

L0DirBouts = [];
AloneUndirBouts = [];
for i = BirdsToUse(:)',
    L0DirBouts(end+1) = DirBoutNos(i,1);
    AloneUndirBouts(end+1) = UndirBoutNos(i,end);
end
disp(['Median # of L0 bouts = ', num2str(median(L0DirBouts)), '; range = ', num2str(min(L0DirBouts)), ' - ', num2str(max(L0DirBouts))]);
disp(['Median # of undir bouts = ', num2str(median(AloneUndirBouts)), '; range = ', num2str(min(AloneUndirBouts)), ' - ', num2str(max(AloneUndirBouts))]);
Harini_INNum_MotifNumDur_L0Dir_UnDir_Diffs(DirBout_Stats(BirdsToUse), DirUndirBout_Stats(BirdsToUse), UndirBout_Stats(BirdsToUse), Distances, BirdNames(BirdsToUse));


% Now to look at effects of distance
% Use only birds that have data for at least 3 distances
BirdsToUse = find((DirBoutNos(:,1) >= MinTrialNo) & (NumDistances(:) >= 2));
disp(['Number of birds with >= ', num2str(MinTrialNo), ' directed song bouts at L0 and atleast 2 other distances = ', num2str(length(BirdsToUse)), ' / ', num2str(size(DirBoutNos, 1)), ' birds']);
disp(['Median # of distances / bird = ', num2str(median(NumDistances(BirdsToUse))), '; range = ', num2str(min(NumDistances(BirdsToUse))), ' - ', num2str(max(NumDistances(BirdsToUse)))]);

NumBoutsPerDistance = [];
for i = BirdsToUse(:)',
    ValidDistances = find(DirBoutNos(i,:) >= MinTrialNo);
    NumBoutsPerDistance = [NumBoutsPerDistance(:); DirBoutNos(i, ValidDistances)'];
end
disp(['Median # of directed song bouts / distance = ', num2str(median(NumBoutsPerDistance)), '; range = ', num2str(min(NumBoutsPerDistance)), ' - ', num2str(max(NumBoutsPerDistance))]);

switch (BirdOption)
    case 'INNum_MotifNumDur'
        Harini_FeatureFits_INNum_MotifNumDur_Time(DirBout_Stats(BirdsToUse), DirUndirBout_Stats(BirdsToUse), UndirBout_Stats(BirdsToUse), Distances, BirdNames(BirdsToUse));
    case 'Amplitude'
        Harini_FeatureFits_Amplitude_Time(DirBout_Stats(BirdsToUse), DirUndirBout_Stats(BirdsToUse), UndirBout_Stats(BirdsToUse), Distances, BirdNames(BirdsToUse));
    case 'FF'
        Harini_FeatureFits_FF_Time(DirBout_Stats(BirdsToUse), DirUndirBout_Stats(BirdsToUse), UndirBout_Stats(BirdsToUse), Distances, BirdNames(BirdsToUse));
end

disp('Finished calculating bout stats');