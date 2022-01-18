function [] = Harini_AnalyzeSyllAmplitudesPatterns(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

AdjustOrNot = 1; % 0 means use unadjusted syll amplitudes, else 1 means use adjusted syll amplitudes
MinSyllNo = 10;
Fs = 44100; % Sampling rate
FontSizeVal = 12;
OutputDir = '/home/raghav/StudentRelated/Harini';

% Load up the raw waveforms for each bird and calculate means, etc.
for i = 1:length(IndividualBirds),
    disp(BirdNames{i});
    clear TempRawWFs;
    TempRawWFs = load(fullfile('/home/raghav/StudentRelated/Harini', [BirdNames{i}, 'RawSyllLogAmpTraces.mat']));
    
    % Now to get syllable means, medians, etc.
    % 
    % First find all sylls with more than the min syll no - defined at the
    % top
    UniqueSylls = unique(char(IndividualBirds(i).AllSyllableData(:,1)));
    
    SyllsToPlot = [];
    for j = 1:length(UniqueSylls),
        if (~isempty(find(IndividualBirds(i).SortedBirdParameters(1).MotifLabels == UniqueSylls(j))))
            SyllsToPlot(end+1) = j;
        end
    end
    
    UniqueSylls = UniqueSylls(SyllsToPlot);
    
    if (isempty(UniqueSylls))
        disp('No motif syllables');
        continue;
    end
    
    SyllIndex = 0;
    for j = UniqueSylls(:)',
        disp(['Syllable ', j]);
        SyllIndex = SyllIndex + 1;
    
        MatchingSylls = find(char(IndividualBirds(i).AllSyllableData(:,1)) == j);
        MatchingSylls = MatchingSylls(find(MatchingSylls <= length(TempRawWFs.TempAdjustedSyllLogAmplitudes)));

        SessionNum = 0;
        Session(i).Syll(SyllIndex).SyllMean = [];
        Session(i).Syll(SyllIndex).SyllMedian = [];
        Session(i).Syll(SyllIndex).SyllMax = [];
        Session(i).Syll(SyllIndex).RecordingDayIndex = [];
        Session(i).Syll(SyllIndex).ConditionIndex = [];
        Session(i).Syll(SyllIndex).Label = j;
        Session(i).Syll(SyllIndex).SyllTime = [];
        Session(i).Syll(SyllIndex).SyllMeanRMS = [];
        
        for k = 1:length(IndividualBirds(i).RecordingDays),
            UniqueConditions = unique(IndividualBirds(i).AllConditionIndices);
            for ConditionIndex = 1:length(UniqueConditions),
                SessionSylls = find((char(IndividualBirds(i).AllSyllableData(:,1)) == j) & (IndividualBirds(i).AllConditionIndices(:) == UniqueConditions(ConditionIndex)) & (IndividualBirds(i).AllRecordingDayIndices(:) == k));
                
                % First remove all the empty amplitudes
                SessionSyllLengths = cellfun(@length, TempRawWFs.TempAdjustedSyllLogAmplitudes(SessionSylls));
                SessionSylls = SessionSylls(find(SessionSyllLengths > 0));
                
                if (~isempty(SessionSylls))
                    if (AdjustOrNot == 1)
                        Session(i).Syll(SyllIndex).SyllMean{end+1} = cellfun(@mean, TempRawWFs.TempAdjustedSyllLogAmplitudes(SessionSylls));
                        Session(i).Syll(SyllIndex).SyllMedian{end+1} = cellfun(@median, TempRawWFs.TempAdjustedSyllLogAmplitudes(SessionSylls));
                        Session(i).Syll(SyllIndex).SyllMax{end+1} = cellfun(@max, TempRawWFs.TempAdjustedSyllLogAmplitudes(SessionSylls));
                    else
                        Session(i).Syll(SyllIndex).SyllMean{end+1} = cellfun(@mean, TempRawWFs.TempRawSyllLogAmplitudes(SessionSylls));
                        Session(i).Syll(SyllIndex).SyllMedian{end+1} = cellfun(@median, TempRawWFs.TempRawSyllLogAmplitudes(SessionSylls));
                        Session(i).Syll(SyllIndex).SyllMax{end+1} = cellfun(@max, TempRawWFs.TempRawSyllLogAmplitudes(SessionSylls));
                    end
                    Session(i).Syll(SyllIndex).SyllMeanRMS{end+1} = IndividualBirds(i).AllSyllableLogAmplitudeRMS(SessionSylls);
                    Session(i).Syll(SyllIndex).RecordingDayIndex{end+1} = k;
                    Session(i).Syll(SyllIndex).ConditionIndex{end+1} = UniqueConditions(ConditionIndex);
                    Session(i).Syll(SyllIndex).SyllTime{end+1} = IndividualBirds(i).FileTime(IndividualBirds(i).AllSyllableData(SessionSylls, 2)) + IndividualBirds(i).AllSyllableData(SessionSylls,4)/(1000*3600);
                end
            end
        end
    end
end

% % First I'm going to plot the within session variability and the global variability.
% % To check for significance, I will randomise the trials across sessions
% % and then measure the ratio of within-session variability and global
% % variability.
% 
% MeasurementTypes = {'Mean' 'Median' 'Max'};
% for Measure = 1:length(MeasurementTypes),
%     SessionIQR{Measure} = [];
%     GlobalIQR{Measure} = [];
%     RandomRatioCIs{Measure} = [];
%     Significance{Measure} = [];
%     ActualRatio{Measure} = [];
%     
%     for i = 1:length(Session)
%         for j = 1:length(Session(i).Syll),
%             GlobalIQR{Measure}(end+1) = iqr(eval(['cell2mat(Session(', num2str(i), ').Syll(', num2str(j), ').Syll', MeasurementTypes{Measure}, ')']));
%             SessionIQR{Measure}(end+1) = mean(cellfun(@iqr, eval(['Session(', num2str(i), ').Syll(', num2str(j), ').Syll', MeasurementTypes{Measure}])));
%             
%             % Now for each session and syllable, also calculate the 95%
%             % confidence intervals for the ratio of session variability to
%             % global variability
%             for Rep = 1:10000,
%                 RandomAmplitudes = cell2mat(eval(['Session(', num2str(i), ').Syll(', num2str(j), ').Syll', MeasurementTypes{Measure}]));
%                 RandomAmplitudes = RandomAmplitudes(randperm(length(RandomAmplitudes)));
%                 RandomAmplitudes = mat2cell(RandomAmplitudes, 1, cellfun(@length, (eval(['Session(', num2str(i), ').Syll(', num2str(j), ').Syll', MeasurementTypes{Measure}]))));
%                 RandomGlobalIQR(Rep) = iqr(cell2mat(RandomAmplitudes));
%                 RandomSessionIQR(Rep) = mean(cellfun(@iqr, RandomAmplitudes));
%             end
%             RandomRatioCIs{Measure}(end+1,:) = prctile(RandomGlobalIQR./RandomSessionIQR, [2.5 97.5]);
%             clear RandomAmplitudes RandomGlobalIQR RandomSessionIQR;
%             ActualRatio{Measure}(end+1) = GlobalIQR{Measure}(end)/SessionIQR{Measure}(end);
%             if ((ActualRatio{Measure}(end) < RandomRatioCIs{Measure}(end,1)) || (ActualRatio{Measure}(end) > RandomRatioCIs{Measure}(end,2))) 
%                 Significance{Measure}(end+1) = 1;
%             else
%                 Significance{Measure}(end+1) = 0;
%             end
%         end
%     end
% end
% 
% % Now to plot the figures. Will plot session vs. global iqr and will make
% % the circles filled if it is significantly different
% figure;
% set(gcf, 'Color', 'w');
% set(gcf, 'Position', [749 268 1050 400]);
% p = panel();
% p.pack('h', {1/3 1/3 1/3});
%     
% for Measure = 1:length(MeasurementTypes),
%     p(Measure).select();
%     hold on;
%     for i = 1:length(SessionIQR{Measure}),
%         if (Significance{Measure}(i) == 1)
%             plot(SessionIQR{Measure}(i), GlobalIQR{Measure}(i), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
%         else
%             plot(SessionIQR{Measure}(i), GlobalIQR{Measure}(i), 'ko', 'MarkerSize', 6);
%         end
%     end
%     axis tight;
%     Temp = axis;
%     Temp = [0.98*min(Temp) 1.02*max(Temp) 0.98*min(Temp) 1.02*max(Temp)];
%     axis(Temp);
%     plot(Temp(1:2), Temp(1:2), 'k--');
%     if (Measure == 2)
%         xlabel('Mean session IQR');
%     end
%     if (Measure == 1)
%         ylabel('Global IQR');
%     end
%     title(MeasurementTypes{Measure});
% end
% p.marginleft = 20;
% p.margintop = 15;
% p.fontsize = 16;
% p.de.margin = 20;
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print(fullfile(OutputDir, ['AllBirds_SessionIQRvsGlobalIQR.png']), '-dpng');
  
% Now for each bird, I should plot the means and medians from each session.
% Separate plots for each syllable in each bird. Individual figures for
% each bird. Calculate means + 95% CIs or medians + 95% CIs for all of the
% data pooled together. Potentially do an ANOVA across all sessions too
% with a post-hoc analysis
for i = 1:length(Session),
    figure;
    p = panel();
    p.pack(ones(1,length(Session(i).Syll))/length(Session(i).Syll));
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [233 72 1000 300*length(Session(i).Syll)]);
    for j = 1:length(Session(i).Syll),
        p(j).select();
        % Now first find the times for each of the sessions and then plot
        % it out.
        clear RecordingDay_Condition_String RecordingTime SyllMeans SyllMedians SyllIQR SyllSEM;
        UniqueConditions = unique(cell2mat(Session(i).Syll(j).ConditionIndex)),
        for k = 1:length(Session(i).Syll(j).SyllMean),
            ConditionIndex = find(UniqueConditions == Session(i).Syll(j).ConditionIndex{k});
            RecordingDay_Condition_String{k} = [IndividualBirds(i).RecordingDays{Session(i).Syll(j).RecordingDayIndex{k}}, ': ', IndividualBirds(i).Conditions{ConditionIndex}]; 
            RecordingTime(k) = Session(i).Syll(j).SyllTime{k}(1) + (Session(i).Syll(j).RecordingDayIndex{k} - 1)*24;
        end
        [SortedVals, SortedIndices] = sort(RecordingTime);
        AllAmps = [];
        AllGroups = [];
        for k = SortedIndices(:)',
            AllAmps = [AllAmps(:); Session(i).Syll(j).SyllMean{k}(:)];
            AllGroups = [AllGroups(:); ones(size(Session(i).Syll(j).SyllMean{k}(:)))*k];
        end
        boxplot(AllAmps, AllGroups);
        Temp = axis;
        % Now plot the overall median and the overall 95% CIs for the whole
        % data
        hold on;
        plot(Temp(1:2), ones(1,2)*median(AllAmps), 'r--');
        plot(Temp(1:2), ones(1,2)*prctile(AllAmps, 25), 'k--');
        plot(Temp(1:2), ones(1,2)*prctile(AllAmps, 75), 'k--');
        
        if (j == 1)
            title([BirdNames{i}, ': Syll ', Session(i).Syll(j).Label]);
        end
        set(gca, 'XTick', 1:1:length(SortedIndices), 'XTickLabel', RecordingDay_Condition_String(SortedIndices), 'XTickLabelRotation', 30);
        ylabel('Log Amplitude (dB)');
    end
    p.fontsize = 14;
    p.de.margin = 30;
    p.marginleft = 20;
    p.marginbottom = 25;
    set(gcf, 'PaperPositionMode', 'auto');
    print(fullfile(OutputDir, [BirdNames{i}, '.SessionWiseAmplitudes_MeanAronovFee.Boxplot.png']), '-dpng');
end

% Now for each bird, I should plot the means and medians from each session.
% Separate plots for each syllable in each bird. Individual figures for
% each bird. Calculate means + 95% CIs or medians + 95% CIs for all of the
% data pooled together. Potentially do an ANOVA across all sessions too
% with a post-hoc analysis
for i = 1:length(Session),
    figure;
    p = panel();
    p.pack(ones(1,length(Session(i).Syll))/length(Session(i).Syll));
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [233 72 1000 300*length(Session(i).Syll)]);
    for j = 1:length(Session(i).Syll),
        p(j).select();
        % Now first find the times for each of the sessions and then plot
        % it out.
        clear RecordingDay_Condition_String RecordingTime SyllMeans SyllMedians SyllIQR SyllSEM;
        UniqueConditions = unique(cell2mat(Session(i).Syll(j).ConditionIndex)),
        for k = 1:length(Session(i).Syll(j).SyllMean),
            ConditionIndex = find(UniqueConditions == Session(i).Syll(j).ConditionIndex{k});
            RecordingDay_Condition_String{k} = [IndividualBirds(i).RecordingDays{Session(i).Syll(j).RecordingDayIndex{k}}, ': ', IndividualBirds(i).Conditions{ConditionIndex}]; 
            RecordingTime(k) = Session(i).Syll(j).SyllTime{k}(1) + (Session(i).Syll(j).RecordingDayIndex{k} - 1)*24;
        end
        [SortedVals, SortedIndices] = sort(RecordingTime);
        AllAmps = [];
        AllGroups = [];
        for k = SortedIndices(:)',
            AllAmps = [AllAmps(:); Session(i).Syll(j).SyllMeanRMS{k}(:)];
            AllGroups = [AllGroups(:); ones(size(Session(i).Syll(j).SyllMeanRMS{k}(:)))*k];
        end
        boxplot(AllAmps, AllGroups);
        Temp = axis;
        % Now plot the overall median and the overall 95% CIs for the whole
        % data
        hold on;
        plot(Temp(1:2), ones(1,2)*median(AllAmps), 'r--');
        plot(Temp(1:2), ones(1,2)*prctile(AllAmps, 25), 'k--');
        plot(Temp(1:2), ones(1,2)*prctile(AllAmps, 75), 'k--');
        
        if (j == 1)
            title([BirdNames{i}, ': Syll ', Session(i).Syll(j).Label]);
        end
        set(gca, 'XTick', 1:1:length(SortedIndices), 'XTickLabel', RecordingDay_Condition_String(SortedIndices), 'XTickLabelRotation', 30);
        ylabel('RMS Amplitude');
    end
    p.fontsize = 14;
    p.de.margin = 30;
    p.marginleft = 20;
    p.marginbottom = 25;
    set(gcf, 'PaperPositionMode', 'auto');
    print(fullfile(OutputDir, [BirdNames{i}, '.SessionWiseAmplitudes.MeanRMS.Boxplot.png']), '-dpng');
end

disp('Finished');

