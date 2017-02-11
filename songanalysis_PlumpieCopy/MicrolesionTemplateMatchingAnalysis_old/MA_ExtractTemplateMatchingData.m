function [Parameters] = MA_ExtractTemplateMatchingData(Parameters, OutputDir)


FileSep = filesep;
%============= Load results of syllable template matching =================
% Next for the 100% shuffled song comparisons

TemplateMatchOutputDir = fullfile(OutputDir, 'ShuffledSongComparisons');
[SyllableTemplateDir, SyllableTemplateFileName, SyllableTemplateExt] = fileparts(Parameters.SyllableTemplateFileName);

TemplateMatchOutputDir = fullfile(TemplateMatchOutputDir, [SyllableTemplateFileName, SyllableTemplateExt, '.TemplateMatchResults']);
disp('Loading results of shuffled song template matching ...');

NumberofFiles = min(50, length(Parameters.PreUnDirSongFileNames{1}));
RandomFiles = Parameters.PreUnDirSongFileNames{1}(1:NumberofFiles);

for SyllTemp = 1:length(Parameters.SyllableTemplate.SyllableTemplates),
    disp(['   Syllable ', Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label]);
    
    LabelCellArray = cellstr(char(ones(length(RandomFiles), 1)*double(['Syll_', Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label])));
    OutputDirCellArray = cellstr(char(ones(length(RandomFiles), 1)*double([TemplateMatchOutputDir, FileSep, 'Syll_', Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label])));
    
    Parameters.SyllableShuffledSongMatches{SyllTemp} = cellfun(@MA_LoadShuffledTemplateMatchResultsFile_WithBoutLens, RandomFiles, OutputDirCellArray, LabelCellArray, 'UniformOutput', 0);
    Parameters.SyllableShuffledSongMatches{SyllTemp} = sort(cell2mat(Parameters.SyllableShuffledSongMatches{SyllTemp}));
    ShuffledSongMatchesThreshold{SyllTemp} = Parameters.SyllableShuffledSongMatches{SyllTemp}(round(0.99 * length(Parameters.SyllableShuffledSongMatches{SyllTemp})));
end

ShufflePercentages = [25 50 75];

ShuffleIndex = 0;
for Shuffle = ShufflePercentages,

    ShuffleIndex = ShuffleIndex + 1;
    
    TemplateMatchOutputDir = fullfile(OutputDir, 'PartShuffledSongComparisons');

    TemplateMatchOutputDir = fullfile(TemplateMatchOutputDir, [SyllableTemplateFileName, SyllableTemplateExt, '.Shuffle.', num2str(Shuffle), 'percent.TemplateMatchResults']);

    disp(['Loading results of ', num2str(Shuffle), ' % shuffled song template matching ...']);

    NumberofFiles = min(50, length(Parameters.PreUnDirSongFileNames{1}));
    RandomFiles = Parameters.PreUnDirSongFileNames{1}(1:NumberofFiles);

    for SyllTemp = 1:length(Parameters.SyllableTemplate.SyllableTemplates),
        
        disp(['   Syllable ', Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label]);
    
        LabelCellArray = cellstr(char(ones(length(RandomFiles), 1)*double(['Syll_', Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label])));
        OutputDirCellArray = cellstr(char(ones(length(RandomFiles), 1)*double([TemplateMatchOutputDir, FileSep, 'Syll_', Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label])));

        Parameters.SyllablePartShuffledSongMatches{ShuffleIndex}{SyllTemp} = cellfun(@MA_LoadShuffledTemplateMatchResultsFile_WithBoutLens, RandomFiles, OutputDirCellArray, LabelCellArray, 'UniformOutput', 0);
        Parameters.SyllablePartShuffledSongMatches{ShuffleIndex}{SyllTemp} = sort(cell2mat(Parameters.SyllablePartShuffledSongMatches{ShuffleIndex}{SyllTemp}));
        PartShuffledSongMatchesThreshold{ShuffleIndex}{SyllTemp} = Parameters.SyllablePartShuffledSongMatches{ShuffleIndex}{SyllTemp}(round(0.95 * length(Parameters.SyllablePartShuffledSongMatches{ShuffleIndex}{SyllTemp})));
    end
end


% Loading results for pre song with threshold as 0.1th percentile of the
% shuffled song matches values (max value from each file)

TemplateMatchOutputDir = fullfile(OutputDir, [SyllableTemplateFileName, SyllableTemplateExt, '.TemplateMatchResults']);

disp('Loading results of song template matching ...');

for SyllTemp = 1:length(Parameters.SyllableTemplate.SyllableTemplates),
    disp(['   Syllable ', Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label]);
    
    Threshold = ShuffledSongMatchesThreshold{SyllTemp};
    
    disp('Loading results of song template matching ...');
    % Load results for all days with zero threshold - can then modify threshold
    % later
    % First for pre song
    for i = 1:Parameters.NoPreDays,
        % First for dir song
        disp(['   Pre Day #', num2str(i), ' - directed song ...']); 
        
        LabelCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double(['Syll_', Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label])));
        OutputDirCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double([TemplateMatchOutputDir, FileSep, 'Syll_', Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label])));
        ThresholdCellArray = num2cell(ones(length(Parameters.PreDirSongFileNames{i}), 1) * Threshold);
        SongFileNoCellArray = num2cell([1:1:length(Parameters.PreDirSongFileNames{i})]');
        % using cellfun so that i iterate over each element of the cell array.
        % To use cellfun, all of the other inputs also have to be in the form
        % of cell arrays of the same length - so the previous three lines
        % convert file type, data dir and output dir - common parameters for
        % all of the files into cell arrays
        Parameters.SyllablePreDirResults{i}{SyllTemp} = cellfun(@MA_LoadTemplateMatchResultsFile_WithBoutLens, Parameters.PreDirSongFileNames{i}, OutputDirCellArray, ThresholdCellArray, LabelCellArray, SongFileNoCellArray, 'UniformOutput', 0);
        Parameters.SyllablePreDirBoutLengths{i}{SyllTemp} = cellfun(@MA_LoadBoutLength_WithBoutLens, Parameters.PreDirSongFileNames{i}, OutputDirCellArray, LabelCellArray, SongFileNoCellArray, 'UniformOutput', 0);
        
        % next undirected songs
        disp(['   Pre Day #', num2str(i), ' - undirected song ...']); 

        LabelCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double(['Syll_', Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label])));
        OutputDirCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double([TemplateMatchOutputDir, FileSep, 'Syll_', Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label])));
        ThresholdCellArray = num2cell(ones(length(Parameters.PreUnDirSongFileNames{i}), 1) * Threshold);
        SongFileNoCellArray = num2cell([1:1:length(Parameters.PreUnDirSongFileNames{i})]');
        % using cellfun so that i iterate over each element of the cell array.
        % To use cellfun, all of the other inputs also have to be in the form
        % of cell arrays of the same length - so the previous three lines
        % convert file type, data dir and output dir - common parameters for
        % all of the files into cell arrays
        Parameters.SyllablePreUnDirResults{i}{SyllTemp} = cellfun(@MA_LoadTemplateMatchResultsFile_WithBoutLens, Parameters.PreUnDirSongFileNames{i}, OutputDirCellArray, ThresholdCellArray, LabelCellArray, SongFileNoCellArray, 'UniformOutput', 0);
        Parameters.SyllablePreUnDirBoutLengths{i}{SyllTemp} = cellfun(@MA_LoadBoutLength_WithBoutLens, Parameters.PreUnDirSongFileNames{i}, OutputDirCellArray, LabelCellArray, SongFileNoCellArray, 'UniformOutput', 0);
    end

    % Next for post song
    for i = 1:Parameters.NoPostDays,
        % First for dir song
        disp(['   Post Day #', num2str(i), ' - directed song ...']); 
        LabelCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double(['Syll_', Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label])));
        OutputDirCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double([TemplateMatchOutputDir, FileSep, 'Syll_', Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label])));
        ThresholdCellArray = num2cell(ones(length(Parameters.PostDirSongFileNames{i}), 1) * Threshold);
        SongFileNoCellArray = num2cell([1:1:length(Parameters.PostDirSongFileNames{i})]');
        % using cellfun so that i iterate over each element of the cell array.
        % To use cellfun, all of the other inputs also have to be in the form
        % of cell arrays of the same length - so the previous three lines
        % convert file type, data dir and output dir - common parameters for
        % all of the files into cell arrays
        Parameters.SyllablePostDirResults{i}{SyllTemp} = cellfun(@MA_LoadTemplateMatchResultsFile_WithBoutLens, Parameters.PostDirSongFileNames{i}, OutputDirCellArray, ThresholdCellArray, LabelCellArray, SongFileNoCellArray, 'UniformOutput', 0);
        Parameters.SyllablePostDirBoutLengths{i}{SyllTemp} = cellfun(@MA_LoadBoutLength_WithBoutLens, Parameters.PostDirSongFileNames{i}, OutputDirCellArray, LabelCellArray, SongFileNoCellArray, 'UniformOutput', 0);

        % next undirected songs
        disp(['   Post Day #', num2str(i), ' - undirected song ...']); 

        LabelCellArray = cellstr(char(ones(length(Parameters.PostUnDirSongFileNames{i}), 1)*double(['Syll_', Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label])));
        OutputDirCellArray = cellstr(char(ones(length(Parameters.PostUnDirSongFileNames{i}), 1)*double([TemplateMatchOutputDir, FileSep, 'Syll_', Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label])));
        ThresholdCellArray = num2cell(ones(length(Parameters.PostUnDirSongFileNames{i}), 1) * Threshold);
        SongFileNoCellArray = num2cell([1:1:length(Parameters.PostUnDirSongFileNames{i})]');
        % using cellfun so that i iterate over each element of the cell array.
        % To use cellfun, all of the other inputs also have to be in the form
        % of cell arrays of the same length - so the previous three lines
        % convert file type, data dir and output dir - common parameters for
        % all of the files into cell arrays
        Parameters.SyllablePostUnDirResults{i}{SyllTemp} = cellfun(@MA_LoadTemplateMatchResultsFile_WithBoutLens, Parameters.PostUnDirSongFileNames{i}, OutputDirCellArray, ThresholdCellArray, LabelCellArray, SongFileNoCellArray, 'UniformOutput', 0);
        Parameters.SyllablePostUnDirBoutLengths{i}{SyllTemp} = cellfun(@MA_LoadBoutLength_WithBoutLens, Parameters.PostUnDirSongFileNames{i}, OutputDirCellArray, LabelCellArray, SongFileNoCellArray, 'UniformOutput', 0);
    end
end

TemplateMatchOutputDir = fullfile(OutputDir, [SyllableTemplateFileName, SyllableTemplateExt, '.TemplateMatchResults']);
save(fullfile(TemplateMatchOutputDir, [SyllableTemplateFileName, SyllableTemplateExt, '.TemplateMatchResults.SavedData.mat']), 'Parameters', 'PartShuffledSongMatchesThreshold', 'ShuffledSongMatchesThreshold');
%==========================================================================