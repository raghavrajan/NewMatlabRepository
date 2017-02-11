function [] = MA_PlotTemplateMatchValues(ParameterFile, OutputDir, ThresholdContext, ThresholdDetermination, Syllables, AllowedErrorPercentage)

Context{1} = 'UnDir';
Context{2} = 'Dir';

FileSep = filesep;
Colours = 'rgbcmk';

%====== Load and extract parameters =======================================
disp('Extracting parameters ...');
Parameters = MA_ParseParametersFile(ParameterFile);
%==========================================================================

%======Now extract all the song file names=================================
disp('Extracting song file names ...');
% First for the pre-treatment days
for i = 1:Parameters.NoPreDays,
    disp(['Extracting song files for pre-treatment day #', num2str(i), ' ...']);
    Parameters.PreDirSongFileNames{i} = MA_ExtractSongFileNames(Parameters.PreDirSongFileList{i});
    Parameters.PreUnDirSongFileNames{i} = MA_ExtractSongFileNames(Parameters.PreUnDirSongFileList{i});
end

% Next for the post-treatment days
for i = 1:Parameters.NoPostDays,
    disp(['Extracting song files for post-treatment day #', num2str(i), ' ...']);    
    Parameters.PostDirSongFileNames{i} = MA_ExtractSongFileNames(Parameters.PostDirSongFileList{i});
    Parameters.PostUnDirSongFileNames{i} = MA_ExtractSongFileNames(Parameters.PostUnDirSongFileList{i});
end
%==========================================================================

%====================== Loading templates =================================
disp('Loading motif template ...');
Parameters.SyllableTemplate = load(Parameters.SyllableTemplateFileName);

%==========================================================================

%==========================================================================
[SyllableTemplateDir, SyllableTemplateFileName, SyllableTemplateExt] = fileparts(Parameters.SyllableTemplateFileName);
TemplateMatchOutputDir = fullfile(OutputDir, [SyllableTemplateFileName, SyllableTemplateExt, '.TemplateMatchResults']);
SavedDataFile = fullfile(TemplateMatchOutputDir, [SyllableTemplateFileName, SyllableTemplateExt, '.TemplateMatchResults.SavedData.mat']);

if (~exist(SavedDataFile, 'file'))
    [Parameters] = MA_ExtractTemplateMatchingData(Parameters, OutputDir);
end
load(SavedDataFile);

%======Now segment files (Aronov Fee style) ===============================
% Segment each file separately, minimum interval = 7ms, minimum duration =
% 7ms and the thresholds are determined separately for each file similar to
% Aronov and Fee (J.Neurosci) paper.

NoteFileDir = '/home/raghav/MicrolesionAnalysisResults/';
FileSep = filesep;

disp('Segmenting files into note files ...');

% First for pre days
for i = 1:Parameters.NoPreDays,
    disp(['   Pre Day #', num2str(i), ' - directed song ...']); 
    % First directed songs
    FileTypeCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double(Parameters.FileType)));
    RawDataDirCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double(Parameters.PreDataDir{i})));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double([NoteFileDir, Parameters.BirdName, '.NoteFiles', FileSep])));
    
    % using cellfun so that i iterate over each element of the cell array.
    % To use cellfun, all of the other inputs also have to be in the form
    % of cell arrays of the same length - so the previous three lines
    % convert file type, data dir and output dir - common parameters for
    % all of the files into cell arrays
    
    [Parameters.PreDirOnsets{i}, Parameters.PreDirOffsets{i}, Parameters.PreDirSyllDurs{i}, Parameters.PreDirGapDurs{i}, Parameters.PreDirThresholds{i}, Parameters.PreDirLens{i}, Parameters.PreDirFs{i}] = cellfun(@MA_SegmentFiles, Parameters.PreDirSongFileNames{i}, FileTypeCellArray, RawDataDirCellArray, OutputDirCellArray, 'UniformOutput', 0);

    % next undirected songs
    disp(['   Pre Day #', num2str(i), ' - undirected song ...']); 
    FileTypeCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double(Parameters.FileType)));
    RawDataDirCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double(Parameters.PreDataDir{i})));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double([NoteFileDir, Parameters.BirdName, '.NoteFiles', FileSep])));
    [Parameters.PreUnDirOnsets{i}, Parameters.PreUnDirOffsets{i}, Parameters.PreUnDirSyllDurs{i}, Parameters.PreUnDirGapDurs{i}, Parameters.PreUnDirThresholds{i}, Parameters.PreUnDirLens{i}, Parameters.PreUnDirFs{i}] = cellfun(@MA_SegmentFiles, Parameters.PreUnDirSongFileNames{i}, FileTypeCellArray, RawDataDirCellArray, OutputDirCellArray, 'UniformOutput', 0);
end

% Next for post days
for i = 1:Parameters.NoPostDays,
    % First directed songs
    disp(['   Post Day #', num2str(i), ' - directed song ...']); 
    FileTypeCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double(Parameters.FileType)));
    RawDataDirCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double(Parameters.PostDataDir{i})));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double([NoteFileDir, Parameters.BirdName, '.NoteFiles', FileSep])));
    
    % using cellfun so that i iterate over each element of the cell array.
    % To use cellfun, all of the other inputs also have to be in the form
    % of cell arrays of the same length - so the previous three lines
    % convert file type, data dir and output dir - common parameters for
    % all of the files into cell arrays
    
    [Parameters.PostDirOnsets{i}, Parameters.PostDirOffsets{i}, Parameters.PostDirSyllDurs{i}, Parameters.PostDirGapDurs{i}, Parameters.PostDirThresholds{i}, Parameters.PostDirLens{i}, Parameters.PostDirFs{i}] = cellfun(@MA_SegmentFiles, Parameters.PostDirSongFileNames{i}, FileTypeCellArray, RawDataDirCellArray, OutputDirCellArray, 'UniformOutput', 0);

    % next undirected songs
    disp(['   Post Day #', num2str(i), ' - undirected song ...']); 
    FileTypeCellArray = cellstr(char(ones(length(Parameters.PostUnDirSongFileNames{i}), 1)*double(Parameters.FileType)));
    RawDataDirCellArray = cellstr(char(ones(length(Parameters.PostUnDirSongFileNames{i}), 1)*double(Parameters.PostDataDir{i})));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PostUnDirSongFileNames{i}), 1)*double([NoteFileDir, Parameters.BirdName, '.NoteFiles', FileSep])));
    [Parameters.PostUnDirOnsets{i}, Parameters.PostUnDirOffsets{i}, Parameters.PostUnDirSyllDurs{i}, Parameters.PostUnDirGapDurs{i}, Parameters.PostUnDirThresholds{i}, Parameters.PostUnDirLens{i}, Parameters.PostUnDirFs{i}] = cellfun(@MA_SegmentFiles, Parameters.PostUnDirSongFileNames{i}, FileTypeCellArray, RawDataDirCellArray, OutputDirCellArray, 'UniformOutput', 0);
end
%==========================================================================

%====================== Loading templates =================================
disp('Loading motif template ...');
Parameters.SyllableTemplate = load(Parameters.SyllableTemplateFileName);

%==========================================================================

% Edges = 0:0.1:10;
% 
% %=========== Plotting the results =========================================
% for SyllTemp = 1:length(Parameters.SyllableTemplate.SyllableTemplates),
%     figure;
%     AllMatchValues = [];
%     Legend{SyllTemp} = [];
%     for i = 1:Parameters.NoPreDays,
%        Temp = [];
%        for j = 1:length(Parameters.SyllablePreDirResults{i}{SyllTemp}),
%            if (~isempty(Parameters.SyllablePreDirResults{i}{SyllTemp}{j}))
%                Temp = [Temp(:); Parameters.SyllablePreDirResults{i}{SyllTemp}{j}(:,1)];
%            end
%        end
%        subplot(2,1,1);
%        hold on;
%        plot(Edges, histc(Temp, Edges)/length(Temp), Colours(i));
%        AllMatchValues = [AllMatchValues; Temp(:)];
%        Legend{SyllTemp}{end+1} = ['Pre Day #', num2str(i)];
%        
%        Temp = [];
%        for j = 1:length(Parameters.SyllablePreUnDirResults{i}{SyllTemp}),
%            if (~isempty(Parameters.SyllablePreUnDirResults{i}{SyllTemp}{j}))
%                Temp = [Temp(:); Parameters.SyllablePreUnDirResults{i}{SyllTemp}{j}(:,1)];
%            end
%        end
%        
%        subplot(2,1,2);
%        hold on;
%        plot(Edges, histc(Temp, Edges)/length(Temp), Colours(i));
%        AllMatchValues = [AllMatchValues; Temp(:)];
%     end
%     
%     for i = 1:Parameters.NoPostDays,
%        Temp = [];
%        for j = 1:length(Parameters.SyllablePostDirResults{i}{SyllTemp}),
%            if (~isempty(Parameters.SyllablePostDirResults{i}{SyllTemp}{j}))
%                Temp = [Temp(:); Parameters.SyllablePostDirResults{i}{SyllTemp}{j}(:,1)];
%            end
%        end
%        
%        subplot(2,1,1);
%        hold on;
%        plot(Edges, histc(Temp, Edges)/length(Temp), Colours(i + Parameters.NoPreDays));
%        AllMatchValues = [AllMatchValues; Temp(:)];
%        Legend{SyllTemp}{end+1} = ['Post Day #', num2str(i)];
%        
%        Temp = [];
%        for j = 1:length(Parameters.SyllablePostUnDirResults{i}{SyllTemp}),
%            if (~isempty(Parameters.SyllablePostUnDirResults{i}{SyllTemp}{j}))
%                Temp = [Temp(:); Parameters.SyllablePostUnDirResults{i}{SyllTemp}{j}(:,1)];
%            end
%        end
%        AllMatchValues = [AllMatchValues; Temp(:)];
%        
%        subplot(2,1,2);
%        hold on;
%        plot(Edges, histc(Temp, Edges)/length(Temp), Colours(i + Parameters.NoPreDays));
%     end
%     subplot(2,1,1);
%     axis tight;
%     Temp1 = axis;
%     subplot(2,1,2);
%     Temp2 = axis;
%     
%     Temp = [Edges(1) max(AllMatchValues)*1.05 0 1.02*max(Temp1(4), Temp2(4))];
%     
%     % Now plot thresholds
%     for j = 1:2,
%         subplot(2,1,j);
%         axis(Temp);
%         set(gca, 'FontSize', 14);
%         
%         plot(ShuffledSongMatchesThreshold{SyllTemp}*ones(1,2), Temp(3:4), 'k--');
%         if (j == 1)
%             Legend{SyllTemp}{end+1} = '100% shuffled song';
%         end
%         
%         for i = 1:length(PartShuffledSongMatchesThreshold),
%             plot(PartShuffledSongMatchesThreshold{i}{SyllTemp}*ones(1,2), Temp(3:4), [Colours(i), '--']);
%             if (j == 1)
%                 Legend{SyllTemp}{end+1} = [num2str(i*25), '% shuffled song'];
%             end
%         end
%         
%         switch j
%             case 1
%                 title(['Syll ', Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label, ' : Directed song'], 'FontSize', 16);
%                 
%             case 2
%                 title(['Syll ', Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label, ' : Undirected song'], 'FontSize', 16);
%                 xlabel('Template Match Value', 'FontSize', 16);
%         end
%         
%         ylabel('Fraction of matches', 'FontSize', 16);
%         
%     end
%     set(gcf, 'Color', 'w');
%     set(gcf, 'Position', [150 32 800 650]);
%     subplot(2,1,1);
%     legend(Legend{SyllTemp});
% end

%=========== Plotting spectrograms ========================================
% Like Glaze and Troyer 2013, I'm going to pull out all the template match
% values between the 25% shuffled song matches threshold and 75% shuffled
% song matches threshold. I will then divide these up into 0.1 size bins
% and make sure there are atleast 20 examples in each of the bins. Then I
% will plot out the spectrorams for these 20 and ask the user to decide on
% how many are correctly classified and how many are not. Once this is
% repeated for the whole set, I can generalize this to the overall data and
% choose a threshold that minimizes false positive and false negative
% rates.
% This process will be done for each of the days on one of the two contexts
% - I can do this for directed song and apply the same thresholds to
% undirected song. I can do vice versa too.


if (ThresholdDetermination == 1)
    
    % Save thresholds, false positives, false negatives into file
    [SyllableTemplateDir, SyllableTemplateFileName, SyllableTemplateExt] = fileparts(Parameters.SyllableTemplateFileName);
    TemplateMatchOutputDir = fullfile(OutputDir, [SyllableTemplateFileName, SyllableTemplateExt, '.TemplateMatchResults']);
    SavedDataFile = fullfile(TemplateMatchOutputDir, [SyllableTemplateFileName, SyllableTemplateExt, '.TemplateMatchResults.OptimalThresholds.', Context{1}, '.', Context{2}, '.SavedData.mat']);
    
    if (exist(SavedDataFile, 'file'))
        load(SavedDataFile);
    
        for SyllTemp = 1:length(Parameters.SyllableTemplate.SyllableTemplates),

            TemplateLen{SyllTemp} = size(Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(23).MotifTemplate, 2) * Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(23).FFTWinOverlap/1000;

        end
    else
        for SyllTemp = 1:length(Parameters.SyllableTemplate.SyllableTemplates),

            TemplateLen{SyllTemp} = size(Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(23).MotifTemplate, 2) * Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(23).FFTWinOverlap/1000;

            for i = 1:Parameters.NoPreDays,
                for k = 1:length(Context),
                    Temp = [];
                    for j = 1:length(eval(['Parameters.SyllablePre', Context{k}, 'Results{i}{SyllTemp}'])),
                        if (~isempty(eval(['Parameters.SyllablePre', Context{k}, 'Results{i}{SyllTemp}{j}'])))
                            Temp = [Temp; eval(['Parameters.SyllablePre', Context{k}, 'Results{i}{SyllTemp}{j}'])];
                        end
                    end
                    if (~isempty(strfind('Dir', Context{k})))
                        Combine = 0;
                    else
                        Combine = 1;
                    end
                    [CorrectMatches_Threshold{k}{SyllTemp}{i}] = MA_ManualSpectClassificationCheck(Temp, PartShuffledSongMatchesThreshold{3}{SyllTemp}, PartShuffledSongMatchesThreshold{1}{SyllTemp}*1.25, Parameters, i, SyllTemp, TemplateLen{SyllTemp}, 'Pre', Context{k}, Combine);
                end
            end
        
            for i = 1:Parameters.NoPostDays,
                for k = 1:length(Context),
                    Temp = [];
                    for j = 1:length(eval(['Parameters.SyllablePost', Context{k}, 'Results{i}{SyllTemp}'])),
                        if (~isempty(eval(['Parameters.SyllablePost', Context{k}, 'Results{i}{SyllTemp}{j}'])))
                            Temp = [Temp; eval(['Parameters.SyllablePost', Context{k}, 'Results{i}{SyllTemp}{j}'])];
                        end
                    end
                    if (~isempty(strfind('Dir', Context{k})))
                        Combine = 0;
                    else
                        Combine = 1;
                    end
                    if (~isempty(Temp))
                        [CorrectMatches_Threshold{k}{SyllTemp}{i + Parameters.NoPreDays}] = MA_ManualSpectClassificationCheck(Temp, PartShuffledSongMatchesThreshold{3}{SyllTemp}*1.25, PartShuffledSongMatchesThreshold{1}{SyllTemp}, Parameters, i, SyllTemp, TemplateLen{SyllTemp}, 'Post', Context{k}, Combine);
                    end
                end                    
            end
        end

        % Save thresholds, false positives, false negatives into file
        [SyllableTemplateDir, SyllableTemplateFileName, SyllableTemplateExt] = fileparts(Parameters.SyllableTemplateFileName);
        TemplateMatchOutputDir = fullfile(OutputDir, [SyllableTemplateFileName, SyllableTemplateExt, '.TemplateMatchResults']);
        SavedDataFile = fullfile(TemplateMatchOutputDir, [SyllableTemplateFileName, SyllableTemplateExt, '.TemplateMatchResults.OptimalThresholds.', Context{1}, '.', Context{2}, '.SavedData.mat']);
        save(SavedDataFile, 'CorrectMatches_Threshold');
    end    

    % Now estimate the false positive and false negative percentages
    for SyllTemp = 1:length(Parameters.SyllableTemplate.SyllableTemplates),
        for k = 1:length(Context),
            for i = 1:Parameters.NoPreDays,
                [Thresholds{k}{SyllTemp}{i}, ErrorRate{k}{SyllTemp}{i}, NonLinearErrorRate{k}{SyllTemp}{i}, OptimalThreshold_BothContext{k}{SyllTemp}{i}] = MA_EstimateFalsePositiveFalseNegativePercentage(CorrectMatches_Threshold{k}{SyllTemp}{i}, AllowedErrorPercentage/100);
            end
        end
    
        for i = 1:Parameters.NoPostDays,
            if (length(CorrectMatches_Threshold{k}{SyllTemp}) >= (i + Parameters.NoPreDays))
                [Thresholds{k}{SyllTemp}{i + Parameters.NoPreDays}, ErrorRate{k}{SyllTemp}{i + Parameters.NoPreDays}, NonLinearErrorRate{k}{SyllTemp}{i + Parameters.NoPreDays}, OptimalThreshold_BothContext{k}{SyllTemp}{i + Parameters.NoPreDays}] = MA_EstimateFalsePositiveFalseNegativePercentage(CorrectMatches_Threshold{k}{SyllTemp}{i + Parameters.NoPreDays}, AllowedErrorPercentage/100);
            end
        end
    end
    
    % Now take optimal thresholds only from one condition based on input
    % argument
    for SyllTemp = 1:length(Parameters.SyllableTemplate.SyllableTemplates),
        for i = 1:Parameters.NoPreDays,
            if (~isempty(strfind(ThresholdContext, Context{1})))
                OptimalThreshold{SyllTemp}{i} = OptimalThreshold_BothContext{1}{SyllTemp}{i};
            else
                OptimalThreshold{SyllTemp}{i} = OptimalThreshold_BothContext{2}{SyllTemp}{i};
            end
        end
        
        for i = 1:Parameters.NoPostDays,
            if (~isempty(strfind(Context{1}, ThresholdContext)))
                OptimalThreshold{SyllTemp}{i + Parameters.NoPreDays} = OptimalThreshold_BothContext{1}{SyllTemp}{i + Parameters.NoPreDays};
            else
                OptimalThreshold{SyllTemp}{i + Parameters.NoPreDays} = OptimalThreshold_BothContext{2}{SyllTemp}{i + Parameters.NoPreDays};
            end
        end
    end
    
    % Now plot the mean match value for each syllable for all days
    
    for SyllTemp = 1:length(Parameters.SyllableTemplate.SyllableTemplates),
        SyllLabel(SyllTemp) = Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label;
    end
        
    for SyllTemp = 1:length(Parameters.SyllableTemplate.SyllableTemplates),
        for i = 1:Parameters.NoPreDays,
            [PreDir{SyllTemp}{i}] = MA_GetMeanandConsistencyForSyllables(Parameters.SyllablePreDirResults{i}{SyllTemp}, OptimalThreshold{SyllTemp}{i}, Parameters.SyllablePreDirBoutLengths{i}{SyllTemp});
            [PreUnDir{SyllTemp}{i}] = MA_GetMeanandConsistencyForSyllables(Parameters.SyllablePreUnDirResults{i}{SyllTemp}, OptimalThreshold{SyllTemp}{i}, Parameters.SyllablePreUnDirBoutLengths{i}{SyllTemp});
        end
        for i = 1:Parameters.NoPostDays,
            if (sum(cellfun(@length,Parameters.SyllablePostDirResults{i}{SyllTemp})) ~= 0)
                [PostDir{SyllTemp}{i}] = MA_GetMeanandConsistencyForSyllables(Parameters.SyllablePostDirResults{i}{SyllTemp}, OptimalThreshold{SyllTemp}{i + Parameters.NoPreDays}, Parameters.SyllablePostDirBoutLengths{i}{SyllTemp});
                [PostUnDir{SyllTemp}{i}] = MA_GetMeanandConsistencyForSyllables(Parameters.SyllablePostUnDirResults{i}{SyllTemp}, OptimalThreshold{SyllTemp}{i + Parameters.NoPreDays}, Parameters.SyllablePostUnDirBoutLengths{i}{SyllTemp});
            end
        end
    end
    
    
    for SyllTemp = 1:length(Parameters.SyllableTemplate.SyllableTemplates),
        if (~isempty(find(Syllables == SyllLabel(SyllTemp))))
            for i = 1:Parameters.NoPreDays,
                MeanMatchValue{SyllTemp}(1,i) = PreDir{SyllTemp}{i}.MeanMatchValue(1);
                STDMatchValue{SyllTemp}(1,i) = PreDir{SyllTemp}{i}.MeanMatchValue(3);
                
                MeanMatchValue{SyllTemp}(2,i) = PreUnDir{SyllTemp}{i}.MeanMatchValue(1);
                STDMatchValue{SyllTemp}(2,i) = PreUnDir{SyllTemp}{i}.MeanMatchValue(3);
                
                MeanConsistency{SyllTemp}(1,i) = PreDir{SyllTemp}{i}.Consistency(1);
                STDConsistency{SyllTemp}(1,i) = PreDir{SyllTemp}{i}.Consistency(3);
                
                MeanConsistency{SyllTemp}(2,i) = PreUnDir{SyllTemp}{i}.Consistency(1);
                STDConsistency{SyllTemp}(2,i) = PreUnDir{SyllTemp}{i}.Consistency(3);
                
                MeanIntervalBetweenMatches{SyllTemp}(1,i) = PreDir{SyllTemp}{i}.MeanIntervalBetweenMatches(1);
                MeanIntervalBetweenMatches{SyllTemp}(1,i) = PreDir{SyllTemp}{i}.MeanIntervalBetweenMatches(3);
                
                MeanIntervalBetweenMatches{SyllTemp}(2,i) = PreUnDir{SyllTemp}{i}.MeanIntervalBetweenMatches(1);
                MeanIntervalBetweenMatches{SyllTemp}(2,i) = PreUnDir{SyllTemp}{i}.MeanIntervalBetweenMatches(3);
            end
            
            for i = 1:Parameters.NoPostDays,
                MeanMatchValue{SyllTemp}(1,i + Parameters.NoPreDays) = PostDir{SyllTemp}{i}.MeanMatchValue(1);
                STDMatchValue{SyllTemp}(1,i + Parameters.NoPreDays) = PostDir{SyllTemp}{i}.MeanMatchValue(3);
                
                MeanMatchValue{SyllTemp}(2,i + Parameters.NoPreDays) = PostUnDir{SyllTemp}{i}.MeanMatchValue(1);
                STDMatchValue{SyllTemp}(2,i + Parameters.NoPreDays) = PostUnDir{SyllTemp}{i}.MeanMatchValue(3);
                
                MeanConsistency{SyllTemp}(1,i + Parameters.NoPreDays) = PostDir{SyllTemp}{i}.Consistency(1);
                STDConsistency{SyllTemp}(1,i + Parameters.NoPreDays) = PostDir{SyllTemp}{i}.Consistency(3);
                
                MeanConsistency{SyllTemp}(2,i + Parameters.NoPreDays) = PostUnDir{SyllTemp}{i}.Consistency(1);
                STDConsistency{SyllTemp}(2,i + Parameters.NoPreDays) = PostUnDir{SyllTemp}{i}.Consistency(3);
                
                if (~isempty(PostDir{SyllTemp}{i}.MeanIntervalBetweenMatches))
                    MeanIntervalBetweenMatches{SyllTemp}(1,i + Parameters.NoPreDays) = PostDir{SyllTemp}{i}.MeanIntervalBetweenMatches(1);
                    MeanIntervalBetweenMatches{SyllTemp}(1,i + Parameters.NoPreDays) = PostDir{SyllTemp}{i}.MeanIntervalBetweenMatches(3);
                else
                    MeanIntervalBetweenMatches{SyllTemp}(1,i + Parameters.NoPreDays) = NaN;
                    MeanIntervalBetweenMatches{SyllTemp}(1,i + Parameters.NoPreDays) = NaN;
                end
                if (~isempty(PostUnDir{SyllTemp}{i}.MeanIntervalBetweenMatches))
                    MeanIntervalBetweenMatches{SyllTemp}(2,i + Parameters.NoPreDays) = PostUnDir{SyllTemp}{i}.MeanIntervalBetweenMatches(1);
                    MeanIntervalBetweenMatches{SyllTemp}(2,i + Parameters.NoPreDays) = PostUnDir{SyllTemp}{i}.MeanIntervalBetweenMatches(3);
                else
                    MeanIntervalBetweenMatches{SyllTemp}(2,i + Parameters.NoPreDays) = NaN;
                    MeanIntervalBetweenMatches{SyllTemp}(2,i + Parameters.NoPreDays) = NaN;
                end
            end
        end
    end
    
    % Now to look at sequences
    
    for i = 1:Parameters.NoPreDays,
        Index = 1;
        for SongFile = 1:length(Parameters.SyllablePreDirBoutLengths{1}{1})
            for BoutNo = 1:size(Parameters.SyllablePreDirBoutLengths{1}{1}{SongFile},1),
                PreDirSequences{i}{Index} = [];
                for SyllTemp = 1:length(Parameters.SyllableTemplate.SyllableTemplates),
                    if (~isempty(find(Syllables == SyllLabel(SyllTemp))))
                        Indices = find((PreDir{SyllTemp}{i}.MatchVals(:,3) == SongFile) & (PreDir{SyllTemp}{i}.MatchVals(:,4) == BoutNo));
                        PreDirSequences{i}{Index}((end+1):(end+length(Indices)),:) = [ones(size(Indices(:)))*SyllTemp PreDir{SyllTemp}{i}.MatchVals(Indices,2)];
                    end
                end
                [SortedVals, SortedIndices] = sort(PreDirSequences{i}{Index}(:,2));
                PreDirSequences{i}{Index} = PreDirSequences{i}{Index}(SortedIndices,:);
                Index = Index + 1;
            end
        end
        
        Index = 1;
        for SongFile = 1:length(Parameters.SyllablePreUnDirBoutLengths{1}{1})
            for BoutNo = 1:size(Parameters.SyllablePreUnDirBoutLengths{1}{1}{SongFile},1),
                PreUnDirSequences{i}{Index} = [];
                for SyllTemp = 1:length(Parameters.SyllableTemplate.SyllableTemplates),
                    if (~isempty(find(Syllables == SyllLabel(SyllTemp))))
                        Indices = find((PreUnDir{SyllTemp}{i}.MatchVals(:,3) == SongFile) & (PreUnDir{SyllTemp}{i}.MatchVals(:,4) == BoutNo));
                        PreUnDirSequences{i}{Index}((end+1):(end+length(Indices)),:) = [ones(size(Indices(:)))*SyllTemp PreUnDir{SyllTemp}{i}.MatchVals(Indices,2)];
                    end
                end
                [SortedVals, SortedIndices] = sort(PreUnDirSequences{i}{Index}(:,2));
                PreUnDirSequences{i}{Index} = PreUnDirSequences{i}{Index}(SortedIndices,:);
                Index = Index + 1;
            end
        end
    end

    for i = 1:Parameters.NoPostDays,
        Index = 1;
        for SongFile = 1:length(Parameters.SyllablePostDirBoutLengths{1}{1})
            for BoutNo = 1:size(Parameters.SyllablePostDirBoutLengths{1}{1}{SongFile},1),
                PostDirSequences{i}{Index} = [];
                for SyllTemp = 1:length(Parameters.SyllableTemplate.SyllableTemplates),
                    if (~isempty(find(Syllables == SyllLabel(SyllTemp))))
                        Indices = find((PostDir{SyllTemp}{i}.MatchVals(:,3) == SongFile) & (PostDir{SyllTemp}{i}.MatchVals(:,4) == BoutNo));
                        PostDirSequences{i}{Index}((end+1):(end+length(Indices)),:) = [ones(size(Indices(:)))*SyllTemp PostDir{SyllTemp}{i}.MatchVals(Indices,2)];
                    end
                end
                [SortedVals, SortedIndices] = sort(PostDirSequences{i}{Index}(:,2));
                PostDirSequences{i}{Index} = PostDirSequences{i}{Index}(SortedIndices,:);
                Index = Index + 1;
            end
        end
        
        Index = 1;
        for SongFile = 1:length(Parameters.SyllablePostUnDirBoutLengths{1}{1})
            for BoutNo = 1:size(Parameters.SyllablePostUnDirBoutLengths{1}{1}{SongFile},1),
                PostUnDirSequences{i}{Index} = [];
                for SyllTemp = 1:length(Parameters.SyllableTemplate.SyllableTemplates),
                    if (~isempty(find(Syllables == SyllLabel(SyllTemp))))
                        if (~isempty(PostUnDir{SyllTemp}{i}.MatchVals))
                            Indices = find((PostUnDir{SyllTemp}{i}.MatchVals(:,3) == SongFile) & (PostUnDir{SyllTemp}{i}.MatchVals(:,4) == BoutNo));
                            if (~isempty(Indices))
                                PostUnDirSequences{i}{Index}((end+1):(end+length(Indices)),:) = [ones(size(Indices(:)))*SyllTemp PostUnDir{SyllTemp}{i}.MatchVals(Indices,2)];
                            end
                        end
                    end
                end
                if (~isempty(PostUnDirSequences{i}{Index}))
                    [SortedVals, SortedIndices] = sort(PostUnDirSequences{i}{Index}(:,2));
                    PostUnDirSequences{i}{Index} = PostUnDirSequences{i}{Index}(SortedIndices,:);                
                end
                Index = Index + 1;
            end
        end
    end
end
disp('Finished');