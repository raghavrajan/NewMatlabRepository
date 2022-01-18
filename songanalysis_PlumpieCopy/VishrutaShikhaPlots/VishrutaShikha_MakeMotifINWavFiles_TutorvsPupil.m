function [] = VishrutaShikha_MakeMotifINWavFiles_TutorvsPupil(BirdParameters)

% Initialise random number generator to same state
rng('default');
OutputDirectory = '/home/raghav/VishrutaTutorPupilMotifINWavefiles';
if (~exist(OutputDirectory, 'dir'))
    mkdir(OutputDirectory);
end
if (~exist(fullfile(OutputDirectory, 'INs')))
    mkdir(fullfile(OutputDirectory, 'INs'));
end

for i = 1:length(BirdParameters),
    % First find valid song bouts
    disp(['Bird #', num2str(i)]);
    SongBouts = find((BirdParameters(i).Bouts(:,7) == 1) & (BirdParameters(i).Bouts(:,8) >= 0) & (BirdParameters(i).Bouts(:,9) >= 1));
    
    % Now choose 20 of these bouts at random and make wav files out of the
    % first motif from each of these bouts and make wav files of 1 each of
    % the different types of INs from these bouts.
    
    % Pre and post padding is the silence period before the motif/IN and
    % after the motif/IN
    
    PrePadding = 10; % in ms
    PostPadding = 10; % in ms
    
    RandBoutIndices = randperm(length(SongBouts));
    RandBouts = SongBouts(RandBoutIndices(1:min(length(RandBoutIndices), 20)));
    
    fprintf('\n');
    for j = 1:length(RandBouts),
        fprintf('%d > ', j);
        
        % First get Bout labels
        BoutIndices = find(BirdParameters(i).SyllableListBoutNum == SongBouts(j));
        BoutLabels = char(BirdParameters(i).SyllableData(BoutIndices,1));
        
        % First plot for motifs
        for k = 1:length(BirdParameters(i).CommonMotifs),
            Motifs = strfind(BoutLabels(:)', BirdParameters(i).CommonMotifs{k});
            if (~isempty(Motifs))
                MotifOnsetFile = BirdParameters(i).SyllableData(Motifs(1) + (BoutIndices(1) - 1), 2);
                MotifOffsetFile = BirdParameters(i).SyllableData(Motifs(1) + (length(BirdParameters(i).CommonMotifs{k}) - 1) + (BoutIndices(1) - 1), 3);
                if (MotifOnsetFile == MotifOffsetFile)
                    MotifOnsetTime = BirdParameters(i).SyllableData(Motifs(1) + (BoutIndices(1) - 1), 4) + PrePadding;
                    MotifOffsetTime = BirdParameters(i).SyllableData(Motifs(1) + (length(BirdParameters(i).CommonMotifs{k}) - 1) + (BoutIndices(1) - 1), 5) + PostPadding;
                    
                    [RawSong, Fs] = GetData(BirdParameters(i).DataDirectory, BirdParameters(i).SongFileNames{MotifOnsetFile}, BirdParameters(i).FileType, 0);
                    Motif = [];
                    Flag = 1;
                    try
                        Motif = RawSong(round(Fs * MotifOnsetTime/1000):round(Fs * MotifOffsetTime/1000));
                    catch
                        Flag = 0;
                    end
                    if (Flag == 1)
                        if (Fs ~= 44100)
                            Motif = resample(Motif,44100,Fs);
                        end
                        if (BirdParameters(i).Tutor == 1)
                            OutputFileName = fullfile(OutputDirectory, [BirdParameters(i).BirdName, '.', BirdParameters(i).UndirSongFileList, '.Nest#', num2str(BirdParameters(i).NestNo), '.Tutor.Motif#', num2str(k), '.Motif#', num2str(j), '.wav']);
                        else
                            OutputFileName = fullfile(OutputDirectory, [BirdParameters(i).BirdName, '.', BirdParameters(i).UndirSongFileList, '.Nest#', num2str(BirdParameters(i).NestNo), '.Pupil.Motif#', num2str(k), '.Motif#', num2str(j), '.wav']);
                        end
                        audiowrite(OutputFileName, Motif, Fs);
                        PlotSpectrogram(OutputDirectory, OutputFileName, 'wav', 'hot');
                        axis([0 1 300 8000]);
                        set(gcf, 'Color', 'w');
                        set(gca, 'Visible', 'off');
                        set(gca, 'Position', [0.05 0.05 0.9 0.9]);
                        set(gcf, 'PaperPositionMode', 'auto');
                        print([OutputFileName, '.eps'], '-depsc2', '-r300');
                    end
                end
            end
        end
        
        % Next plot for INs
        for k = 1:length(BirdParameters(i).INLabels),
            INs = strfind(BoutLabels(:)', BirdParameters(i).INLabels(k));
            if (~isempty(INs))
                INOnsetFile = BirdParameters(i).SyllableData(INs(round(length(INs)/2)) + (BoutIndices(1) - 1), 2);
                INOffsetFile = BirdParameters(i).SyllableData(INs(round(length(INs)/2)) + (BoutIndices(1) - 1), 3);
                if (INOnsetFile == INOffsetFile)
                    INOnsetTime = BirdParameters(i).SyllableData(INs(round(length(INs)/2)) + (BoutIndices(1) - 1), 4) + PrePadding;
                    INOffsetTime = BirdParameters(i).SyllableData(INs(round(length(INs)/2)) + (BoutIndices(1) - 1), 5) + PostPadding;
                    
                    [RawSong, Fs] = GetData(BirdParameters(i).DataDirectory, BirdParameters(i).SongFileNames{INOnsetFile}, BirdParameters(i).FileType, 0);
                    IN = [];
                    Flag = 1;
                    try
                        IN = RawSong(round(Fs * INOnsetTime/1000):round(Fs * INOffsetTime/1000));
                    catch
                        Flag = 0;
                    end
                    if (Flag == 1)
                        if (Fs ~= 44100)
                            IN = resample(IN,44100,Fs);
                        end
                        if (BirdParameters(i).Tutor == 1)
                            OutputFileName = fullfile(OutputDirectory, 'INs', [BirdParameters(i).BirdName, '.', BirdParameters(i).UndirSongFileList, '.Nest#', num2str(BirdParameters(i).NestNo), '.Tutor.IN#', num2str(k), '.', BirdParameters(i).INLabels(k), '.IN#', num2str(j), '.wav']);
                        else
                            OutputFileName = fullfile(OutputDirectory, 'INs', [BirdParameters(i).BirdName, '.', BirdParameters(i).UndirSongFileList, '.Nest#', num2str(BirdParameters(i).NestNo), '.Pupil.IN#', num2str(k), '.', BirdParameters(i).INLabels(k), '.IN#', num2str(j), '.wav']);
                        end
                        audiowrite(OutputFileName, IN, Fs);
                        PlotSpectrogram(fullfile(OutputDirectory, 'INs'), OutputFileName, 'wav', 'hot');
                        axis([0 0.125 300 8000]);
                        set(gcf, 'Color', 'w');
                        set(gca, 'Visible', 'off');
                        set(gcf, 'Position', [153 440 100 250]);
                        set(gca, 'Position', [0.05 0.05 0.9 0.9]);
                        set(gcf, 'PaperPositionMode', 'auto');
                        print([OutputFileName, '.eps'], '-depsc2', '-r300');
                        close all;
                    end
                end
            end
        end
    end
end

disp('Finished');

