function [] = Harini_PlotSyllAmplitudes_AcrossSessions(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

MinSyllNo = 10;
Fs = 44100; % Sampling rate
FontSizeVal = 12;
OutputDir = '/home/raghav/StudentRelated/Harini';

% Load up the raw waveforms for each bird and plot them out
for i = 1:length(IndividualBirds),
    disp(BirdNames{i});
    clear TempRawWFs;
    TempRawWFs = load(fullfile('/home/raghav/StudentRelated/Harini', [BirdNames{i}, 'RawSyllLogAmpTraces.mat']));
    
    % Now to plot out the raw amplitudes and means for each session and
    % each day separately. In addition will plot the means across all
    % data as an overlaid trace in all plots. Will have to warp the data
    % for this to the median syllable length
    
    % Will do this for all syllables that have atleast 10 renditions across
    % all of the data
    
    % First find all sylls with more than the min syll no - defined at the
    % top
    UniqueSylls = unique(char(IndividualBirds(i).AllSyllableData(:,1)));
    SyllsToPlot = [];
    for j = 1:length(UniqueSylls),
        NumSylls = length(find(char(IndividualBirds(i).AllSyllableData(:,1)) == UniqueSylls(j)));
        if (NumSylls >= MinSyllNo)
            SyllsToPlot(end+1) = j;
        end
    end
    
    UniqueSylls = UniqueSylls(SyllsToPlot);
    
    % Now find total number of sessions and divide the figure accordingly
    % for each syllable separately
    TotalNumSessions = 0;
    for j = 1:length(IndividualBirds(i).RecordingDays),
        UniqueConditions = unique(IndividualBirds(i).AllConditionIndices);
        for k = UniqueConditions(:)',
            MatchingSessions = find((IndividualBirds(i).AllConditionIndices == k) & (IndividualBirds(i).AllRecordingDayIndices == j));
            if (~isempty(MatchingSessions))
                TotalNumSessions = TotalNumSessions + 1;
            end
        end
    end
    
    % Now make a figure with panels to plot each of the sessions separately
    % Separate figures for each syllable.
    % I will keep column no to a max of 3 and then add the appropriate # of
    % rows
    if ((TotalNumSessions == 0) || isempty(UniqueSylls))
        disp('No sessions OR no unique syllables');
        continue;
    end
    
    NumCols = 3; % # of columns
    
    SyllIndex = 0;
    for j = UniqueSylls(:)',
        disp(['Syllable ', j]);
        SyllIndex = SyllIndex + 1;
        figure;
        p = panel();
        NumRows = ceil(TotalNumSessions/NumCols);
        p.pack(NumRows, NumCols);
    
        % Now for each syll to plot, I should find the median length and warp
        % all renditions to the median. Then calculate overall average and sem
        MatchingSylls = find(char(IndividualBirds(i).AllSyllableData(:,1)) == j);
        MatchingSylls = MatchingSylls(find(MatchingSylls <= length(TempRawWFs.TempAdjustedSyllLogAmplitudes)));
        MedianSyllLength = median(cellfun(@length, TempRawWFs.TempAdjustedSyllLogAmplitudes(MatchingSylls)));
        
        WarpedSyllLogAmplitudes = [];
        % Now warp all sylls to that length
        for k = MatchingSylls(:)',
            if (~isempty(TempRawWFs.TempAdjustedSyllLogAmplitudes{k}))
                WarpedSyllLogAmplitudes(k,:) = spline(1:1:length(TempRawWFs.TempAdjustedSyllLogAmplitudes{k}),TempRawWFs.TempAdjustedSyllLogAmplitudes{k}, linspace(1, length(TempRawWFs.TempAdjustedSyllLogAmplitudes{k}), MedianSyllLength));
            end
        end
        
        if (isempty(WarpedSyllLogAmplitudes))
            disp('All empty log amplitude waveforms');
            continue;
        end
        
        SessionNum = 0;
        Session(i).Syll(SyllIndex).SyllMean = [];
        Session(i).Syll(SyllIndex).RecordingDayIndex = [];
        Session(i).Syll(SyllIndex).ConditionIndex = [];
        
        for k = 1:length(IndividualBirds(i).RecordingDays),
            UniqueConditions = unique(IndividualBirds(i).AllConditionIndices);
            for ConditionIndex = 1:length(UniqueConditions),
                SessionSylls = find((char(IndividualBirds(i).AllSyllableData(:,1)) == j) & (IndividualBirds(i).AllConditionIndices(:) == UniqueConditions(ConditionIndex)) & (IndividualBirds(i).AllRecordingDayIndices(:) == k));
                if (length(SessionSylls) >= 3)
                    Session(i).Syll(SyllIndex).SyllMean{end+1} = mean(WarpedSyllLogAmplitudes(SessionSylls,:),2);
                    Session(i).Syll(SyllIndex).RecordingDayIndex(end+1) = k;
                    Session(i).Syll(SyllIndex).ConditionIndex(end+1) = UniqueConditions(ConditionIndex);
                    
                    SessionNum = SessionNum + 1;
                    RowNum = ceil(SessionNum/NumCols);
                    ColNum = SessionNum - ((RowNum - 1) * NumCols);
                    p(RowNum, ColNum).select();
                    plot((1:1:MedianSyllLength)/Fs, WarpedSyllLogAmplitudes(SessionSylls,:)', 'k', 'Color', [0.7 0.7 0.7]);
                    hold on;
                    patch([(1:1:MedianSyllLength)/Fs fliplr((1:1:MedianSyllLength)/Fs)], [(nanmean(WarpedSyllLogAmplitudes(SessionSylls,:)) + (nanstd(WarpedSyllLogAmplitudes(SessionSylls,:))/sqrt(length(SessionSylls)))) fliplr(nanmean(WarpedSyllLogAmplitudes(SessionSylls,:)) - (nanstd(WarpedSyllLogAmplitudes(SessionSylls,:))/sqrt(length(SessionSylls))))], 'k', 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.2);
                    SessionMeanPlot = plot((1:1:MedianSyllLength)/Fs, nanmean(WarpedSyllLogAmplitudes(SessionSylls,:)), 'b');
                    GlobalMeanPlot = plot((1:1:MedianSyllLength)/Fs, nanmean(WarpedSyllLogAmplitudes(MatchingSylls,:)), 'r');
                    title([BirdNames{i}, ': Syll ', j, ': ', IndividualBirds(i).RecordingDays{k}, ': ', IndividualBirds(i).Conditions{ConditionIndex}]);
                    legend([SessionMeanPlot GlobalMeanPlot], {'Session mean', 'Global mean'}, 'Location', 'southwest');
                    legend('boxoff');
                    axis tight;
                    
                    if (RowNum == NumRows)
                        xlabel('Time (ms)');
                    end
                    if (ColNum == 1)
                        ylabel('Log Amplitude (dB)');
                    end
                end
            end
        end
        p.fontsize = FontSizeVal;
        % p.fontname = 'Times';
        p.margintop = 10;
        p.marginleft = 20;
        p.marginright = 20;

        set(gcf, 'Color', 'w');
        set(gcf, 'Position', [680 37 1000 900]);
        set(gcf, 'PaperPositionMode', 'auto');
        print(fullfile(OutputDir, [BirdNames{i}, ': Syll ', j, '.SessionAmplitudeWFs.png']), '-dpng');
        close all;
    end
end

   

disp('Finished');

