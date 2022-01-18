function [PairWise_Corr, RandomPairWise_Corr] = CalcCorrBetweenINSequences(InterBoutInterval)

InterINInterval = 500; % ms - maximum allowed interval between INs
INTiming_Fs = 1000; % Hz - sampling rate for the 1/0 representation of pre-first motif syllable data
% File name of the file used to get the data from
SongDetailsFile = '/home/raghav/StudentRelated/Harini/Manuscript/JEBManuscript/ScriptsForFigures/ContinuousDataSongAnalysis_BirdDetails_ForHariniDataRelatedAnalysis.csv';
SavedDataDir = '/home/raghav/StudentRelated/Harini/Manuscript/JEBManuscript/SavedData';

%% Get all data
[BirdParameters, Flag] = ProcessSongData_IntoBouts(SongDetailsFile, InterBoutInterval, SavedDataDir);

%% Now to calculate correlation for each bout pre-motif
% Steps to do this
% 1) Find all first motif syllables for each song bout
% 2) Now find the earliest starting time and make a 0,1, amplitude curve
% for the pre-first motif syllable time - 0 refers to silence and 1 refers
% to the presence of a syllable. For this I should only consider the last
% set of syllables with < 500ms gap between them

for i = 1:length(BirdParameters),
    ValidSongBouts = find((BirdParameters(i).Bouts(:,8) > 0) & (BirdParameters(i).Bouts(:,9) > 1) & (BirdParameters(i).Bouts(:,7) == 1));
    MotifSylls = BirdParameters(i).MotifLabels;
    PreFirstMotifSyllTime{i} = ones(size(ValidSongBouts))*NaN;
    for j = 1:length(ValidSongBouts),
        BirdParameters(i).BoutDetails{ValidSongBouts(j)}.FirstMotifSyllIndex = NaN;
        BirdParameters(i).BoutDetails{ValidSongBouts(j)}.FirstMotifSyllOnsetTime = NaN;
        BirdParameters(i).BoutDetails{ValidSongBouts(j)}.PreMotifOnsetTime = NaN;
        
        for k = 1:length(BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Labels),
            if (~isempty(find(MotifSylls == BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Labels(k))))
                BirdParameters(i).BoutDetails{ValidSongBouts(j)}.FirstMotifSyllIndex = k;
                BirdParameters(i).BoutDetails{ValidSongBouts(j)}.FirstMotifSyllOnsetTime = BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Onsets(k);
                InterINGaps = BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Onsets(2:k) - BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Offsets(1:k-1);
                LongGaps = find(InterINGaps >= InterINInterval);
                if (isempty(LongGaps))
                    BirdParameters(i).BoutDetails{ValidSongBouts(j)}.PreMotifOnsetTime = -BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Onsets(k);
                    PreFirstMotifSyllTime{i}(j) = -BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Onsets(k);
                else
                    BirdParameters(i).BoutDetails{ValidSongBouts(j)}.PreMotifOnsetTime = -(BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Onsets(k) - BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Onsets(LongGaps(end) + 1));
                    PreFirstMotifSyllTime{i}(j) = -(BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Onsets(k) - BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Onsets(LongGaps(end) + 1));
                end
                break;
            end
        end
    end
    EarliestTime = prctile(PreFirstMotifSyllTime{i}, 25) - iqr(PreFirstMotifSyllTime{i});
    PreFirstMotifSyll_INRepresentation{i} = zeros(length(ValidSongBouts), ceil(abs(EarliestTime) * INTiming_Fs/1000));
    for j = 1:length(ValidSongBouts),
        for k = BirdParameters(i).BoutDetails{ValidSongBouts(j)}.FirstMotifSyllIndex-1:-1:1,
            if ((BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Onsets(k) - BirdParameters(i).BoutDetails{ValidSongBouts(j)}.FirstMotifSyllOnsetTime) >= EarliestTime)
                StartIndex = round(INTiming_Fs/1000 * abs(BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Onsets(k) - BirdParameters(i).BoutDetails{ValidSongBouts(j)}.FirstMotifSyllOnsetTime - EarliestTime));
                if (StartIndex == 0)
                    StartIndex = 1;
                end
                EndIndex = round(INTiming_Fs/1000 * abs(BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Offsets(k) - BirdParameters(i).BoutDetails{ValidSongBouts(j)}.FirstMotifSyllOnsetTime - EarliestTime));
            else
                if ((BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Offsets(k) - BirdParameters(i).BoutDetails{ValidSongBouts(j)}.FirstMotifSyllOnsetTime) >= EarliestTime)
                    StartIndex = 1;
                    EndIndex = round(INTiming_Fs/1000 * abs(BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Offsets(k) - BirdParameters(i).BoutDetails{ValidSongBouts(j)}.FirstMotifSyllOnsetTime - EarliestTime));
                else
                    StartIndex = NaN;
                    EndIndex = NaN;
                end
            end
            
            if (~isnan(StartIndex))
               PreFirstMotifSyll_INRepresentation{i}(j, StartIndex:EndIndex) = 1;
            end
        end
    end
    % Now calculate pairwise correlations
    PairWise_Corr{i} = [];
    for j = 1:size(PreFirstMotifSyll_INRepresentation{i},1),
        for k = j+1:size(PreFirstMotifSyll_INRepresentation{i},1),
            ST1 = PreFirstMotifSyll_INRepresentation{i}(j,:);
            ST2 = PreFirstMotifSyll_INRepresentation{i}(k,:);
            PairWise_Corr{i}(end+1) = ((ST1 - mean(ST1)) * (ST2 - mean(ST2))')/(norm(ST1 - mean(ST1)) * norm(ST2 - mean(ST2)));
        end
    end
    
    % Now to check random pairwise correlations
    % To do this, I will resort the rows in the pre-first motif syll
    % representation to keep the same duration and same inter-IN interval
    % with just a bit of jitter (+/- 30ms). I will wrap the data around to
    % the beginning if it reaches the end.
    % Jitter is from a uniform distribution between -30ms and +30ms
    Jitter = 50; % in ms
    for j = 1:10000,
        RandomPreFirstMotifSyll_INRepresentation = PreFirstMotifSyll_INRepresentation{i};
        for k = 1:size(RandomPreFirstMotifSyll_INRepresentation,1),
            RandomJitter = (rand * 2 * Jitter) - Jitter;
            Temp = find(RandomPreFirstMotifSyll_INRepresentation(k,:) == 1) + round(RandomJitter * INTiming_Fs/1000);
            Temp(find(Temp > size(RandomPreFirstMotifSyll_INRepresentation,2))) = Temp(find(Temp > size(RandomPreFirstMotifSyll_INRepresentation,2))) - size(RandomPreFirstMotifSyll_INRepresentation,2);
            Temp(find(Temp <= 0)) = Temp(find(Temp <= 0)) + size(RandomPreFirstMotifSyll_INRepresentation,2);

            RandomPreFirstMotifSyll_INRepresentation(k,:) = zeros(1, size(RandomPreFirstMotifSyll_INRepresentation,2));
            RandomPreFirstMotifSyll_INRepresentation(k,Temp) = 1;
        end
        
        % Do the mean subtraction now.
        RandomPreFirstMotifSyll_INRepresentation = RandomPreFirstMotifSyll_INRepresentation - repmat(mean(RandomPreFirstMotifSyll_INRepresentation, 2), 1, size(RandomPreFirstMotifSyll_INRepresentation,2));

        TempRandomPairWise_Corr = [];
        Temp2RandomPairWise_Corr = [];
        
        for k = 1:size(RandomPreFirstMotifSyll_INRepresentation,1),
            for l = k+1:size(RandomPreFirstMotifSyll_INRepresentation,1),
                ST1 = RandomPreFirstMotifSyll_INRepresentation(k,:);
                ST2 = RandomPreFirstMotifSyll_INRepresentation(l,:);
                %TempRandomPairWise_Corr(end+1) = ((ST1 - mean(ST1)) * (ST2 - mean(ST2))')/(norm(ST1 - mean(ST1)) * norm(ST2 - mean(ST2)));
                Temp2RandomPairWise_Corr(end+1) = CalculateCorrBetweenSpikeTrains(ST1, ST2);
            end
        end
        RandomPairWise_Corr{i}(j) = mean(Temp2RandomPairWise_Corr);
    end
        
    
    figure;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [680 128 1150 850]);
    p = panel();
    p.pack({4/5 1/5});
    p(1).select();
    imagesc(linspace(EarliestTime, 0, size(PreFirstMotifSyll_INRepresentation{i},2)),1:1:size(PreFirstMotifSyll_INRepresentation{i},1), PreFirstMotifSyll_INRepresentation{i});
    title([BirdParameters(i).BirdName, ': pairwise IN timing correlation = ', num2str(mean(PairWise_Corr{i})), ' +/- ', num2str(std(PairWise_Corr{i})), ' (95%CI: ', num2str(prctile(RandomPairWise_Corr{i}, 2.5)), ' - ', num2str(prctile(RandomPairWise_Corr{i}, 97.5)), ')']);
    ylabel('Bout Index');
    axis tight;
    p(2).select();
    hold on;
    patch([linspace(EarliestTime, 0, size(PreFirstMotifSyll_INRepresentation{i},2)) fliplr(linspace(EarliestTime, 0, size(PreFirstMotifSyll_INRepresentation{i},2)))], [(mean(PreFirstMotifSyll_INRepresentation{i}) - std(PreFirstMotifSyll_INRepresentation{i})/sqrt(size(PreFirstMotifSyll_INRepresentation{i},1))) fliplr((mean(PreFirstMotifSyll_INRepresentation{i}) + std(PreFirstMotifSyll_INRepresentation{i})/sqrt(size(PreFirstMotifSyll_INRepresentation{i},1))))], 'k', 'FaceColor', 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(linspace(EarliestTime, 0, size(PreFirstMotifSyll_INRepresentation{i},2)),mean(PreFirstMotifSyll_INRepresentation{i}), 'k')
    xlabel('Time relative to motif onset (msec)');
    axis tight;
    ylabel('Sound (1) or Silence (0)');
    p.fontsize = 16;
    p.margintop = 15;
    p.de.margin = 30;
    p.marginbottom = 20;
    p.marginleft = 20;
    set(gcf, 'PaperPositionMode', 'auto');
    OutputDir = '/home/raghav/LabRelated';
    print(fullfile(OutputDir, [BirdParameters(i).UndirSongFileList,'.INTiming_Correlation.png']), '-dpng', '-r300');
end
disp('Finished');