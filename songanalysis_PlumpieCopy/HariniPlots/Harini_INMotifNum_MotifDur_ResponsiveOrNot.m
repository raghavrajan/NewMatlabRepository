function [DirBout_Stats, UndirBout_Stats, DirUndirBout_Stats] = Harini_INMotifNum_MotifDur_ResponsiveOrNot(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

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
        DirFemaleResponseBouts{i}{j} = intersect(IndividualBirds(i).BoutStatistics(strmatch('R', IndividualBirds(i).BoutFemaleResponse, 'exact'), BoutIndexColumnIndex), DirBouts{i}{j});
        DirFemaleNoResponseBouts{i}{j} = intersect(IndividualBirds(i).BoutStatistics(strmatch('NR', IndividualBirds(i).BoutFemaleResponse, 'exact'), BoutIndexColumnIndex), DirBouts{i}{j});
        
        FemaleResponseBouts{i}{j} = intersect(IndividualBirds(i).BoutStatistics(strmatch('R', IndividualBirds(i).BoutFemaleResponse, 'exact'), BoutIndexColumnIndex), Bouts{i}{j});
        FemaleNoResponseBouts{i}{j} = intersect(IndividualBirds(i).BoutStatistics(strmatch('NR', IndividualBirds(i).BoutFemaleResponse, 'exact'), BoutIndexColumnIndex), Bouts{i}{j});
        
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
        if (~isempty(DirFemaleResponseBouts{i}{j}))
            [DirFemaleResponseBout_Stats{i}{j}] = Harini_CalculateBoutINMotifNumbers_Duration_AmpKao(IndividualBirds(i), DirFemaleResponseBouts{i}{j}, MotifLabels, INLabels, CommonMotifs);
        else
            DirFemaleResponseBout_Stats{i}{j} = [];
        end
        if (~isempty(DirFemaleNoResponseBouts{i}{j}))
            [DirFemaleNoResponseBout_Stats{i}{j}] = Harini_CalculateBoutINMotifNumbers_Duration_AmpKao(IndividualBirds(i), DirFemaleNoResponseBouts{i}{j}, MotifLabels, INLabels, CommonMotifs);
        else
            DirFemaleNoResponseBout_Stats{i}{j} = [];
        end
    end
end

% First get number of bouts for all types
for i = 1:length(DirBout_Stats),
    DirBoutNos(i,:) = cellfun(@length, DirBouts{i});
    DirFemaleResponseBoutNos(i,:) = cellfun(@length, DirFemaleResponseBouts{i});
    DirFemaleNoResponseBoutNos(i,:) = cellfun(@length, DirFemaleNoResponseBouts{i});
    
    FemaleResponseBoutNos(i,:) = cellfun(@length, FemaleResponseBouts{i});
    FemaleNoResponseBoutNos(i,:) = cellfun(@length, FemaleNoResponseBouts{i});
end

% I need to find outlier thresholds for both motif duration and amplitude
% based on all of the data
for i = 1:length(IndividualBirds),
    CommonMotifs = IndividualBirds(i).SortedBirdParameters(1).CommonMotifs;
    MotifLabels = IndividualBirds(i).SortedBirdParameters(1).MotifLabels;
    ValidSongBouts = find((IndividualBirds(i).Bouts(:,7) == 1) & (IndividualBirds(i).Bouts(:,8) > 0) & (IndividualBirds(i).Bouts(:,9) > 1));
    AllMotifDur{i} = [];
    for k = 1:length(MotifLabels),
        AllSyllAmp{i}{k} = [];
        AllSyllDur{i}{k} = [];
    end
    for k = ValidSongBouts(:)',
        BoutLabels = char(IndividualBirds(i).AllSyllableData(IndividualBirds(i).Bouts(k,1):IndividualBirds(i).Bouts(k,2), 1));
        BoutOnsets = IndividualBirds(i).AllSyllableData(IndividualBirds(i).Bouts(k,1):IndividualBirds(i).Bouts(k,2), 4);
        BoutOffsets = IndividualBirds(i).AllSyllableData(IndividualBirds(i).Bouts(k,1):IndividualBirds(i).Bouts(k,2), 5);
        
        BoutSAPDur = BoutOffsets - BoutOnsets;
        BoutSAPDur = BoutSAPDur(:);
        
        BoutSAPLogAmp = IndividualBirds(i).AllSyllableLogAmplitudeKao(IndividualBirds(i).Bouts(k,1):IndividualBirds(i).Bouts(k,2));
        BoutSAPLogAmp = BoutSAPLogAmp(:);
        
        CompleteMotifs = strfind(BoutLabels(:)', CommonMotifs{1});
        if (~isempty(CompleteMotifs))
            MotifDurations = BoutOffsets(CompleteMotifs + length(CommonMotifs{1}) - 1) - BoutOnsets(CompleteMotifs);
            AllMotifDur{i} = [AllMotifDur{i}; MotifDurations(:)];
        end
        
        for j = 1:length(MotifLabels),
            MatchingSylls = find(BoutLabels == MotifLabels(j));
            if (~isempty(MatchingSylls))
                LogAmplitude = BoutSAPLogAmp(MatchingSylls(:));
                SyllDuration = BoutSAPDur(MatchingSylls(:));
                
                AllSyllAmp{i}{j} = [AllSyllAmp{i}{j}; LogAmplitude(:)];
                AllSyllDur{i}{j} = [AllSyllDur{i}{j}; SyllDuration(:)];
            end
        end
    end
    
    MotifDurOutlierThresholds(i,:) = [(prctile(AllMotifDur{i}, 25) - 3*iqr(AllMotifDur{i})) (prctile(AllMotifDur{i}, 75) + 3*iqr(AllMotifDur{i}))];
    for k = 1:length(MotifLabels),
        SyllDurOutlierThresholds{i}(k,:) = [(prctile(AllSyllDur{i}{k}, 25) - 3*iqr(AllSyllDur{i}{k})) (prctile(AllSyllDur{i}{k}, 75) + 3*iqr(AllSyllDur{i}{k}))];
        ValidAmplitudes = AllSyllAmp{i}{k}(find((AllSyllDur{i}{k} >= SyllDurOutlierThresholds{i}(k,1)) & (AllSyllDur{i}{k} <= SyllDurOutlierThresholds{i}(k,2))));
        SyllAmpOutlierThresholds{i}(k,:) = [(prctile(ValidAmplitudes, 25) - 3*iqr(ValidAmplitudes)) (prctile(ValidAmplitudes, 75) + 3*iqr(ValidAmplitudes))];
    end
end
    
    
MinTrialNo = 3;
MinMotifNo = 7;
INNum_Comparison = [];
MotifNum_Comparison = [];
MotifDuration_Comparison = [];
SyllAmplitude_Comparison = [];

for i = 1:size(DirFemaleResponseBoutNos,1),
    ValidBirds = find(((DirFemaleResponseBoutNos(i,:) >= MinTrialNo) + (DirFemaleNoResponseBoutNos(i,:) >= MinTrialNo)) == 2);
    if (~isempty(ValidBirds))
        FemaleResponseINNum = 0;
        FemaleNoResponseINNum = 0;
        FemaleResponseMotifNum = 0;
        FemaleNoResponseMotifNum = 0;
        
        for j = ValidBirds(:)'
            FemaleResponseINNum = FemaleResponseINNum + nanmean(DirFemaleResponseBout_Stats{i}{j}.TotalINNumber_500ms);
            FemaleNoResponseINNum = FemaleNoResponseINNum + nanmean(DirFemaleNoResponseBout_Stats{i}{j}.TotalINNumber_500ms);
            
            FemaleResponseMotifNum = FemaleResponseMotifNum + nanmean(DirFemaleResponseBout_Stats{i}{j}.CompleteMotifNumber);
            FemaleNoResponseMotifNum = FemaleNoResponseMotifNum + nanmean(DirFemaleNoResponseBout_Stats{i}{j}.CompleteMotifNumber);
        end

        INNum_Comparison(end+1,:) = [FemaleResponseINNum/length(ValidBirds) FemaleNoResponseINNum/length(ValidBirds)];
        MotifNum_Comparison(end+1,:) = [FemaleResponseMotifNum/length(ValidBirds) FemaleNoResponseMotifNum/length(ValidBirds)];
        
    end
    
    ValidBirds = find(((DirFemaleResponseBoutNos(i,:) > 0) + (DirFemaleNoResponseBoutNos(i,:) > 0)) == 2);
    if (~isempty(ValidBirds))
        FemaleResponseMotifDur = 0;
        FemaleNoResponseMotifDur = 0;
        
        ValidBirdIndex = 0;
        for j = ValidBirds(:)'
            if ((length(DirFemaleResponseBout_Stats{i}{j}.AllMotifDuration) >= MinMotifNo) && (length(DirFemaleNoResponseBout_Stats{i}{j}.AllMotifDuration) >= MinMotifNo))
                ValidBirdIndex = ValidBirdIndex + 1;
                TempFemaleResponseMotifDurs = DirFemaleResponseBout_Stats{i}{j}.AllMotifDuration;
                TempFemaleNoResponseMotifDurs = DirFemaleNoResponseBout_Stats{i}{j}.AllMotifDuration;
                
                TempFemaleResponseMotifDurs = TempFemaleResponseMotifDurs(find((TempFemaleResponseMotifDurs >= MotifDurOutlierThresholds(i,1)) & (TempFemaleResponseMotifDurs <= MotifDurOutlierThresholds(i,2))));
                TempFemaleNoResponseMotifDurs = TempFemaleNoResponseMotifDurs(find((TempFemaleNoResponseMotifDurs >= MotifDurOutlierThresholds(i,1)) & (TempFemaleNoResponseMotifDurs <= MotifDurOutlierThresholds(i,2))));
                FemaleResponseMotifDur = FemaleResponseMotifDur + nanmean(TempFemaleResponseMotifDurs);
                FemaleNoResponseMotifDur = FemaleNoResponseMotifDur + nanmean(TempFemaleNoResponseMotifDurs);
            end
        end
        if (ValidBirdIndex > 0)
            MotifDuration_Comparison(end+1,:) = [FemaleResponseMotifDur/ValidBirdIndex FemaleNoResponseMotifDur/ValidBirdIndex];
        end
        
        % Now for syllable amplitude
        ValidBirdIndex = 0;
                
        for j = ValidBirds(:)',
            FemaleResponseNumSylls = cellfun(@length, DirFemaleResponseBout_Stats{i}{j}.LogAmplitude);
            FemaleNoResponseNumSylls = cellfun(@length, DirFemaleNoResponseBout_Stats{i}{j}.LogAmplitude);
            ValidSylls = find(((FemaleResponseNumSylls >= MinMotifNo) + (FemaleNoResponseNumSylls >= MinMotifNo)) == 2);
            
            if (j == ValidBirds(1))
                clear FemaleResponseSyllAmplitude FemaleNoResponseSyllAmplitude ValidSyllIndex;
                for k = 1:length(ValidSylls),
                    FemaleResponseSyllAmplitude(k) = 0;
                    FemaleNoResponseSyllAmplitude(k) = 0;
                    ValidSyllIndex(k) = 0;
                end
            end
            for k = ValidSylls(:)',
                ValidSyllIndex(k) = ValidSyllIndex(k) + 1;
                TempFemaleResponseAmplitudes = DirFemaleResponseBout_Stats{i}{j}.LogAmplitude{k};
                TempFemaleResponseDurs = DirFemaleResponseBout_Stats{i}{j}.SyllDuration{k}*1000;
                
                % Check for outliers first in syll duration and then in
                % amplitude
                TempFemaleResponseAmplitudes = TempFemaleResponseAmplitudes(find((TempFemaleResponseDurs >= SyllDurOutlierThresholds{i}(k,1)) & (TempFemaleResponseDurs <= SyllDurOutlierThresholds{i}(k,2)))); 
                TempFemaleResponseAmplitudes = TempFemaleResponseAmplitudes(find((TempFemaleResponseAmplitudes >= SyllAmpOutlierThresholds{i}(k,1)) & (TempFemaleResponseAmplitudes <= SyllAmpOutlierThresholds{i}(k,2)))); 
                
                TempFemaleNoResponseAmplitudes = DirFemaleNoResponseBout_Stats{i}{j}.LogAmplitude{k};
                TempFemaleNoResponseDurs = DirFemaleNoResponseBout_Stats{i}{j}.SyllDuration{k}*1000;
                % Check for outliers first in syll duration and then in
                % amplitude
                TempFemaleNoResponseAmplitudes = TempFemaleNoResponseAmplitudes(find((TempFemaleNoResponseDurs >= SyllDurOutlierThresholds{i}(k,1)) & (TempFemaleNoResponseDurs <= SyllDurOutlierThresholds{i}(k,2)))); 
                TempFemaleNoResponseAmplitudes = TempFemaleNoResponseAmplitudes(find((TempFemaleNoResponseAmplitudes >= SyllAmpOutlierThresholds{i}(k,1)) & (TempFemaleNoResponseAmplitudes <= SyllAmpOutlierThresholds{i}(k,2)))); 
                
                FemaleResponseSyllAmplitude(k) = FemaleResponseSyllAmplitude(k) + nanmean(TempFemaleResponseAmplitudes);
                FemaleNoResponseSyllAmplitude(k) = FemaleNoResponseSyllAmplitude(k) + nanmean(TempFemaleNoResponseAmplitudes);
                
                [Hyp, PVal] = ttest2(TempFemaleResponseAmplitudes, TempFemaleNoResponseAmplitudes);
            end
        end
        for k = 1:length(ValidSyllIndex),
            if (ValidSyllIndex(k) > 0)
                SyllAmplitude_Comparison(end+1,:) = [FemaleResponseSyllAmplitude(k)/ValidSyllIndex(k) FemaleNoResponseSyllAmplitude(k)/ValidSyllIndex(k) i k PVal];
            end
        end
    end
end
            
Distances = [0 20 60 110 165 195];
DistanceLabelString = {'0' '20' '60' '110' '165' 'UN'};

% Now plot the percentage of female response song bouts at each distance

PercentageFemaleResponses = 100 * FemaleResponseBoutNos./(FemaleResponseBoutNos + FemaleNoResponseBoutNos);
DistanceMatrix = repmat(Distances(:)', size(PercentageFemaleResponses,1), 1);
NonNaNPercentages = find(~isnan(PercentageFemaleResponses));
[Corr_R, Corr_P] = corr(DistanceMatrix(NonNaNPercentages), PercentageFemaleResponses(NonNaNPercentages));

Symbols = 'sod^v><+hx*p';


switch BirdOption

    case 'INMotifNum_MotifDur'
        p = panel();
        p.pack({1/2 1/2});
        p(1).pack('h', {1/2 1/2});
        p(2).pack('h', {1/3 1/3 1/3});
        p(1,1).select();
        hold on;
        plot(repmat(Distances, size(PercentageFemaleResponses, 1), 1)', PercentageFemaleResponses', 'ko-', 'Color', [0.75 0.75 0.75]);
        for i = 1:size(PercentageFemaleResponses, 2),
            Indices = find(~isnan(PercentageFemaleResponses(:,i)));
            errorbar(Distances(i), mean(PercentageFemaleResponses(Indices,i)), std(PercentageFemaleResponses(Indices,i))/sqrt(length(Indices)), 'ks', 'MarkerSize', 8, 'LineWidth', 2);
        end
        text(177, 0, '//', 'FontSize', 16)
        axis([-5 Distances(end)+5 0 115]);
        set(gca, 'Box', 'off')
        set(gcf, 'Color', 'w')
        xlabel({'Distance from the male (cm)'});
        set(gca, 'XTick', Distances, 'XTickLabel', DistanceLabelString);
        ylabel({'Proportion of bouts'; 'with female responses (%)'});
        text(Distances(round(length(Distances)/2)+1), 115, {['r = ', num2str(Corr_R)]; ['p = ', num2str(Corr_P)]});

        p(2,1).select();
        hold on;
        INNumBar = bar(mean(INNum_Comparison));
        set(INNumBar, 'FaceColor', 'none');
        for i = 1:size(INNum_Comparison, 2),
            errorbar(i, mean(INNum_Comparison(:,i)), std(INNum_Comparison(:,i))/sqrt(size(INNum_Comparison,1)), 'k');
        end
        plot(repmat([1.15 1.85], size(INNum_Comparison, 1), 1)', INNum_Comparison', 'ko-', 'MarkerSize', 6);
        PValue = signrank(INNum_Comparison(:,1), INNum_Comparison(:,2));
        disp(['P value for comparison of IN number between dir bouts with female responses and without female responses = ', num2str(PValue)]);
        axis tight;
        Temp = axis;
        Temp = [0.5 2.5 0 1.02*Temp(4)];
        axis(Temp);
        sigstar({[1 2]}, PValue);
        ylabel('# of INs');
        set(gca, 'XTick', [1 2], 'XTickLabel', {'Female response' 'No female response'}, 'XTickLabelRotation', 45);

        p(2,2).select();
        hold on;
        MotifNumBar = bar(mean(MotifNum_Comparison));
        set(MotifNumBar, 'FaceColor', 'none');
        for i = 1:size(MotifNum_Comparison, 2),
            errorbar(i, mean(MotifNum_Comparison(:,i)), std(MotifNum_Comparison(:,i))/sqrt(size(MotifNum_Comparison,1)), 'k');
        end
        plot(repmat([1.15 1.85], size(MotifNum_Comparison, 1), 1)', MotifNum_Comparison', 'ko-', 'MarkerSize', 6);
        PValue = signrank(MotifNum_Comparison(:,1), MotifNum_Comparison(:,2));
        disp(['P value for comparison of motif number between dir bouts with female responses and without female responses = ', num2str(PValue)]);
        axis tight;
        Temp = axis;
        Temp = [0.5 2.5 0 1.02*Temp(4)];
        axis(Temp);
        sigstar({[1 2]}, PValue);
        ylabel('# of motifs/bout');
        set(gca, 'XTick', [1 2], 'XTickLabel', {'Female response' 'No female response'}, 'XTickLabelRotation', 45);

        p(2,3).select();
        hold on;
        MotifDurationBar = bar(mean(MotifDuration_Comparison));
        set(MotifDurationBar, 'FaceColor', 'none');
        for i = 1:size(MotifDuration_Comparison, 2),
            errorbar(i, mean(MotifDuration_Comparison(:,i)), std(MotifDuration_Comparison(:,i))/sqrt(size(MotifDuration_Comparison,1)), 'k');
        end
        plot(repmat([1.15 1.85], size(MotifDuration_Comparison, 1), 1)', MotifDuration_Comparison', 'ko-', 'MarkerSize', 6);
        PValue = signrank(MotifDuration_Comparison(:,1), MotifDuration_Comparison(:,2));
        disp(['P value for comparison of motif duration between dir bouts with female responses and without female responses = ', num2str(PValue)]);
        axis tight;
        Temp = axis;
        Temp = [0.5 2.5 0.99*Temp(3) 1.02*Temp(4)];
        axis(Temp);
        sigstar({[1 2]}, PValue);
        ylabel('Motif duration (msec)');
        set(gca, 'XTick', [1 2], 'XTickLabel', {'Female response' 'No female response'}, 'XTickLabelRotation', 45);

        p.de.margin = 22;
        p.fontsize = 12;
        % p.fontname = 'Times';
        p.margintop = 15;
        p.marginleft = 25;
        % p.marginright = 20;
        p.marginbottom = 40;
        set(gcf, 'Color', 'w');
        set(gcf, 'PaperPositionMode', 'auto');
        set(gcf, 'Units', 'inches');
        set(gcf, 'Position', [10 1 8.3 7]);

        print(fullfile(FinalFigureDir, ['FemaleResponseDependence.', BirdOption, '.eps']), '-depsc2', '-r600');
        print(fullfile(FinalFigureDir, ['FemaleResponseDependence', BirdOption, '.png']), '-dpng', '-r600');
        
    case 'Amplitude'
        p = panel();
        p.pack({1});
        p(1).select();
        hold on;
        for i = 1:size(SyllAmplitude_Comparison,1),
            if (SyllAmplitude_Comparison(i,end) < 0.05)
                MarkerSymbolFaceColor = 'k';
            else
                MarkerSymbolFaceColor = 'none';
            end
            plot(SyllAmplitude_Comparison(i,3), SyllAmplitude_Comparison(i,1)/SyllAmplitude_Comparison(i,2), 'ko', 'MarkerFaceColor', 'none', 'MarkerSize', 8, 'LineWidth', 2);
        end
        BirdFinalIndices = [8 10 11 12 13 14 5];

        UniqueBirdIndices = unique(SyllAmplitude_Comparison(:,3));
        for i = 1:length(UniqueBirdIndices),
            BirdIndexXLabelString{i} = ['Bird #', num2str(BirdFinalIndices(UniqueBirdIndices(i)))];
        end
        set(gca, 'XTick', UniqueBirdIndices, 'XTickLabel', BirdIndexXLabelString, 'XTickLabelRotation', 45);
        
        ylabel({'Syll Amplitude Ratio '; 'Female response bouts /'; 'No female response bouts'});
        axis tight;
        Temp = axis;
        Temp = [Temp(1)-0.25 Temp(2)+0.25 0.98*Temp(3) 1.02*Temp(4)];
        axis(Temp);
        plot(Temp(1:2), ones(1,2), 'k--');
        
        p.fontsize = 12;
        % p.fontname = 'Times';
%         p.margintop = 15;
        p.marginleft = 30;
%         % p.marginright = 20;
        p.marginbottom = 20;
        set(gcf, 'Color', 'w');
        set(gcf, 'PaperPositionMode', 'auto');
        set(gcf, 'Units', 'inches');
        set(gcf, 'Position', [10 1 3.5 2.8]);

        print(fullfile(FinalFigureDir, ['SyllAmplitude_FemaleResponseDependence.', BirdOption, '.eps']), '-depsc2', '-r600');
        print(fullfile(FinalFigureDir, ['Syll_Amplitude_FemaleResponseDependence', BirdOption, '.png']), '-dpng', '-r600');
        
end
close all;

disp('Finished calculating bout stats');