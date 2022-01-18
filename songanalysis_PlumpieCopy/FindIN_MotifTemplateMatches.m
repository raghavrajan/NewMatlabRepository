function [] = FindIN_MotifTemplateMatches(BirdDetailsTextFile)

OutputDir = '/data/raghav/ForShikha_INMotif_TemplateMatchAnalysis';
% First get details from the CSV text file
disp('Getting header data from CSV file ...');
[HeaderLine, BirdDetails] = LSINA_GetDetailsFromCSVFile(BirdDetailsTextFile);

% Now parse all the lines into the appropriate variables based on the
% header line
disp('Getting data from CSV file ...');
[BirdParameters] = LSINA_ParseHeaderBirdData(HeaderLine, BirdDetails);

% For each of the birds, now construct templates for the INs and then match
% them to the motifs in the file
for i = 1:length(BirdParameters),
    INLabels = BirdParameters(i).INLabels;
    Notes = load([fullfile(BirdParameters(i).DataDirectory, BirdParameters(i).NoteFileDirectory, BirdParameters(i).SongFileName), '.not.mat']);
    MakeAllSyllableTemplatesFromFile(BirdParameters(i).DataDirectory, BirdParameters(i).NoteFileDirectory, BirdParameters(i).SongFileName, BirdParameters(i).FileType, OutputDir, 0, -5:1:5, 0, setdiff(Notes.labels, INLabels));
    Template = load([fullfile(OutputDir, BirdParameters(i).SongFileName), '.SyllTemplates.0.template.mat']);
    
    [RawData, Fs] = GetData(BirdParameters(i).DataDirectory, BirdParameters(i).SongFileName, BirdParameters(i).FileType, 0);
    Time = (1:1:length(RawData))/Fs;
    [P, F, S1, T] = CalculateMultiTaperSpectrogram(RawData, Fs, 8, 8/2, 1.5);
    
    Freq1 = find((F >= 860) & (F <= 8600));
    S = log10(abs(S1(Freq1,:)));
    % Also, make a 75% shuffled version of the above spectrogram
    % Shuffling only within a time bin
    ShufflePercentage = 75;
    for j = 1:size(S,2),
        % First decide which set to shuffle
        NumToShuffle = round(ShufflePercentage * size(S,1)/100);
        Rand_Shuffle_Indices = randperm(size(S,1));
        RowsToShuffle_Indices = Rand_Shuffle_Indices(1:NumToShuffle);
        RowsNotToShuffle_Indices = Rand_Shuffle_Indices(NumToShuffle+1:end);
        RowsToShuffle = [];
        for k = RowsToShuffle_Indices(:)',
            RowsToShuffle = [RowsToShuffle; S(k,j)];
        end
        RowsToShuffle = RowsToShuffle(randperm(length(RowsToShuffle)));
        Index = 1;
        for k = 1:size(S,1),
            if (~isempty(find(RowsNotToShuffle_Indices == k)))
                RandShuffled_S(k,j) = S(k,j);
            else
                RandShuffled_S(k,j) = RowsToShuffle(Index);
                Index = Index + 1;
            end
        end
    end
    
    SyllableMatchFs = 1/(T(2) - T(1));
    
    
    PrePostPadding = 0.01; % in seconds
    for k = 1:length(Template.SyllableTemplates),
        TempSyllTemplate = Template.SyllableTemplates{k}{min(length(Template.SyllableTemplates{k}),3)};
        fprintf('>');
        TempMatch = ASSLTemplateMatch(S, TempSyllTemplate.MotifTemplate);
        
        % Remove the portion that actually contains the template itself
        Is = find(Notes.labels == TempSyllTemplate.MotifTemplate(1).Label);
        I_To_Remove = min(3, length(Is));
        I_To_Remove = Is(I_To_Remove);
        
        I_To_Remove_OnsetTime = Notes.onsets(I_To_Remove)/1000 - PrePostPadding;
        I_To_Remove_OffsetTime = Notes.offsets(I_To_Remove)/1000 + PrePostPadding;
        
        RandomTempMatch = ASSLTemplateMatch(RandShuffled_S, TempSyllTemplate.MotifTemplate);
        TempMatchTime = linspace(T(1), T(length(TempMatch)), length(TempMatch));

        TempMatch(find((TempMatchTime >= I_To_Remove_OnsetTime) & (TempMatchTime <= I_To_Remove_OffsetTime))) = prctile(TempMatch, 10);
        RandomTempMatch(find((TempMatchTime >= I_To_Remove_OnsetTime) & (TempMatchTime <= I_To_Remove_OffsetTime))) = prctile(RandomTempMatch, 10);
        
        % Now average the motif related match values
        Motifs = strfind(Notes.labels, BirdParameters(i).CommonMotifs{1});
        
    
        Index = 1;
        for j = Motifs(:)',
            MotifOnsetTime = Notes.onsets(j)/1000 - PrePostPadding;
            MotifOffsetTime = Notes.offsets(j + length(BirdParameters(i).CommonMotifs{1}) - 1)/1000 + PrePostPadding;
            MotifDuration(Index) = Notes.offsets(j + length(BirdParameters(i).CommonMotifs{1}) - 1)/1000 - Notes.onsets(j)/1000;
            MotifMatchValues{Index} = TempMatch(find((TempMatchTime >= MotifOnsetTime) & (TempMatchTime <= MotifOffsetTime)));
            RandomMotifMatchValues{Index} = RandomTempMatch(find((TempMatchTime >= MotifOnsetTime) & (TempMatchTime <= MotifOffsetTime)));
            Index = Index + 1;
        end
        
        MinLen = min(cellfun(@length, MotifMatchValues));
        clear MatchValues;
        for j = 1:length(MotifMatchValues),
            MatchValues(j,:) = MotifMatchValues{j}(1:MinLen);
        end
        
        figure;
        p = panel();
        p.pack('h', [1/5 4/5]);
        p(1).pack([1/2 1/2]);
        p(2).pack([1/2 1/2]);
        
        p(1,1).select();
        PlotSpectrogramInAxis_SongVar(RawData(round(I_To_Remove_OnsetTime * Fs):round(I_To_Remove_OffsetTime * Fs)), Time(round(I_To_Remove_OnsetTime * Fs):round(I_To_Remove_OffsetTime * Fs)), Fs, gca);
        set(gca, 'YTick', [0:2000:8000], 'YTickLabel', 0:2:8);
        
        p(2,1).select();
        MotifOnsetTime = Notes.onsets(Motifs(1))/1000 - PrePostPadding;
        MotifOffsetTime = Notes.offsets(Motifs(1) + length(BirdParameters(i).CommonMotifs{1}) - 1)/1000 + PrePostPadding;
        PlotSpectrogramInAxis_SongVar(RawData(round(MotifOnsetTime * Fs):round(MotifOffsetTime * Fs)), Time(round(MotifOnsetTime * Fs):round(MotifOffsetTime * Fs)), Fs, gca);
        ylabel('');
        set(gca, 'YColor', 'w');
        title([BirdParameters(i).SongFileName, ': Motif ', BirdParameters(1).CommonMotifs{1}]);
        
        p(2,2).select();
        hold on;
        MatchValuesTime = linspace(Time(round(MotifOnsetTime * Fs)), Time(round(MotifOffsetTime * Fs)), MinLen);
        patch([MatchValuesTime fliplr(MatchValuesTime)], [(mean(MatchValues) + std(MatchValues)/sqrt(size(MatchValues,1))) fliplr(mean(MatchValues) - std(MatchValues)/sqrt(size(MatchValues,1)))], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        plot(MatchValuesTime, mean(MatchValues), 'b', 'LineWidth', 1);
        plot(MatchValuesTime([1 end]), ones(1,2) * (mean(TempMatch) + 3*std(TempMatch)), 'k--', 'LineWidth', 1);
        axis tight;
        title(['Template match to Intro Note ', TempSyllTemplate.MotifTemplate(1).Label]);
        set(gcf, 'Color', 'w');
        set(gcf, 'Position', [300 200 1100 550]);
        set(gcf, 'PaperPositionMode', 'auto');
        xlabel('Time (sec)');
        
        p.de.margin = 25;
        p.margin = 10;
        p.marginleft = 25;
        p.marginbottom = 20;
        p.fontsize = 14;
        print([BirdParameters(i).BirdName, '.TemplateMatch_IN_to_Motif.png'], '-dpng', '-r300'); 
    end
end

disp('Finished');