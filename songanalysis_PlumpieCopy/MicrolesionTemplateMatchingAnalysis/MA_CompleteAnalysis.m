function [Results, Parameters] = MA_CompleteAnalysis(ParameterFile, OtherBirdsSongParametersFile, Multiplier, ShuffledMatchesFigure, NoteFileDir)

%==========================================================================
% Function for analysing the effects of a treatment on song - written in
% the context of the analysis of effects of HVC microlesions.
% Raghav Rajan - 29th November 2013
%==========================================================================

%============== Some common variables =====================================
OutputDir = '/home/raghav/MicrolesionAnalysisResults/';

FileSep = filesep;
if (OutputDir(end) ~= FileSep)
    OutputDir(end+1) = FileSep;
end
%==========================================================================

%====== Load and extract parameters =======================================
disp('Extracting parameters ...');
Parameters = MA_ParseParametersFile(ParameterFile);
if (~exist([OutputDir, Parameters.BirdName, '.NoteFiles'], 'dir'))
    disp('Creating output directory ...');
    mkdir(OutputDir, [Parameters.BirdName, '.NoteFiles']);
end
    
%==========================================================================

%====Now load up the parameters for the other bird songs===================
% This is a file that is provided as input to the program - it contains
% details of directory and song files for other bird songs
Parameters.OtherBirdsSongDetails = MA_ParseParametersFile(OtherBirdsSongParametersFile);
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

% Next for the other birds songs
disp(['Extracting song files for other birds songs ', num2str(i), ' ...']);
for i = 1:Parameters.OtherBirdsSongDetails.NoPreDays,
    Parameters.OtherBirdsSongDetails.PreDirSongFileNames{i} = MA_ExtractSongFileNames(Parameters.OtherBirdsSongDetails.PreDirSongFileList{i});
    Parameters.OtherBirdsSongDetails.PreUnDirSongFileNames{i} = MA_ExtractSongFileNames(Parameters.OtherBirdsSongDetails.PreUnDirSongFileList{i});
    
    for j = 1:length(Parameters.ExcludeBirds),
        ExcludeFiles = find(cellfun(@length, strfind(Parameters.OtherBirdsSongDetails.PreDirSongFileNames{i}, Parameters.ExcludeBirds{j})));
        if (~isempty(ExcludeFiles))
            Parameters.OtherBirdsSongDetails.PreDirSongFileNames{i}(ExcludeFiles) = [];
            Parameters.OtherBirdsSongDetails.PreUnDirSongFileNames{i}(ExcludeFiles) = [];
        end
    end
end


%==========================================================================

%======Now segment files (Aronov Fee style) ===============================
% Segment each file separately, minimum interval = 7ms, minimum duration =
% 7ms and the thresholds are determined separately for each file similar to
% Aronov and Fee (J.Neurosci) paper.

disp('Segmenting files into note files ...');

% First for pre days
for i = 1:Parameters.NoPreDays,
    disp(['   Pre Day #', num2str(i), ' - directed song ...']); 
    % First directed songs
    FileTypeCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double(Parameters.FileType)));
    RawDataDirCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double(Parameters.PreDataDir{i})));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double([OutputDir, Parameters.BirdName, '.NoteFiles', FileSep])));
    
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
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double([OutputDir, Parameters.BirdName, '.NoteFiles', FileSep])));
    [Parameters.PreUnDirOnsets{i}, Parameters.PreUnDirOffsets{i}, Parameters.PreUnDirSyllDurs{i}, Parameters.PreUnDirGapDurs{i}, Parameters.PreUnDirThresholds{i}, Parameters.PreUnDirLens{i}, Parameters.PreUnDirFs{i}] = cellfun(@MA_SegmentFiles, Parameters.PreUnDirSongFileNames{i}, FileTypeCellArray, RawDataDirCellArray, OutputDirCellArray, 'UniformOutput', 0);
end

% Next for post days
for i = 1:Parameters.NoPostDays,
    % First directed songs
    disp(['   Post Day #', num2str(i), ' - directed song ...']); 
    FileTypeCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double(Parameters.FileType)));
    RawDataDirCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double(Parameters.PostDataDir{i})));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double([OutputDir, Parameters.BirdName, '.NoteFiles', FileSep])));
    
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
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PostUnDirSongFileNames{i}), 1)*double([OutputDir, Parameters.BirdName, '.NoteFiles', FileSep])));
    [Parameters.PostUnDirOnsets{i}, Parameters.PostUnDirOffsets{i}, Parameters.PostUnDirSyllDurs{i}, Parameters.PostUnDirGapDurs{i}, Parameters.PostUnDirThresholds{i}, Parameters.PostUnDirLens{i}, Parameters.PostUnDirFs{i}] = cellfun(@MA_SegmentFiles, Parameters.PostUnDirSongFileNames{i}, FileTypeCellArray, RawDataDirCellArray, OutputDirCellArray, 'UniformOutput', 0);
end
%==========================================================================

%======Now load up note files with bout information =======================
% Load up note files which have bouts labelled with a 'B'. This has been
% done to get a better idea of consistency of song production.

disp('Loading up note files ...');

if (~isempty(NoteFileDir))
    % First for pre days
    for i = 1:Parameters.NoPreDays,
        disp(['   Pre Day #', num2str(i), ' - directed song ...']); 
        % First directed songs
        SlashIndex = find(Parameters.PreDataDir{i} == FileSep);
        
        NoteFileDirCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double([NoteFileDir, FileSep, Parameters.PreDataDir{i}(SlashIndex(end)+1:end)])));

        % using cellfun so that i iterate over each element of the cell array.
        % To use cellfun, all of the other inputs also have to be in the form
        % of cell arrays of the same length - so the previous three lines
        % convert file type, data dir and output dir - common parameters for
        % all of the files into cell arrays

        [Parameters.PreDirBoutOnsets{i}, Parameters.PreDirBoutOffsets{i}, Parameters.PreDirBoutLens{i}] = cellfun(@MA_LoadBoutNoteFiles, Parameters.PreDirSongFileNames{i}, NoteFileDirCellArray, 'UniformOutput', 0);

        % next undirected songs
        disp(['   Pre Day #', num2str(i), ' - undirected song ...']); 
        NoteFileDirCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double([NoteFileDir, FileSep, Parameters.PreDataDir{i}(SlashIndex(end)+1:end)])));
        [Parameters.PreUnDirBoutOnsets{i}, Parameters.PreUnDirBoutOffsets{i}, Parameters.PreUnDirBoutLens{i}] = cellfun(@MA_LoadBoutNoteFiles, Parameters.PreUnDirSongFileNames{i}, NoteFileDirCellArray, 'UniformOutput', 0);
    end

    % Next for post days
    for i = 1:Parameters.NoPostDays,
        % First directed songs
        disp(['   Post Day #', num2str(i), ' - directed song ...']); 
        
        SlashIndex = find(Parameters.PostDataDir{i} == FileSep);
        NoteFileDirCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double([NoteFileDir, FileSep, Parameters.PostDataDir{i}(SlashIndex(end)+1:end)])));

        % using cellfun so that i iterate over each element of the cell array.
        % To use cellfun, all of the other inputs also have to be in the form
        % of cell arrays of the same length - so the previous three lines
        % convert file type, data dir and output dir - common parameters for
        % all of the files into cell arrays

        [Parameters.PostDirBoutOnsets{i}, Parameters.PostDirBoutOffsets{i}, Parameters.PostDirBoutLens{i}] = cellfun(@MA_LoadBoutNoteFiles, Parameters.PostDirSongFileNames{i}, NoteFileDirCellArray, 'UniformOutput', 0);

        % next undirected songs
        disp(['   Post Day #', num2str(i), ' - undirected song ...']); 
        NoteFileDirCellArray = cellstr(char(ones(length(Parameters.PostUnDirSongFileNames{i}), 1)*double([NoteFileDir, FileSep, Parameters.PostDataDir{i}(SlashIndex(end)+1:end)])));
        [Parameters.PostUnDirBoutOnsets{i}, Parameters.PostUnDirBoutOffsets{i}, Parameters.PostUnDirBoutLens{i}] = cellfun(@MA_LoadBoutNoteFiles, Parameters.PostUnDirSongFileNames{i}, NoteFileDirCellArray, 'UniformOutput', 0);
    end
end
%==========================================================================

%======Now plot syllable and gap duration histograms ======================
disp('Plotting syllable and gap duration histograms ...');
Colours = 'rgbcmy';
Symbols = 'o*';

BinSize = 5;
MaxSyllDur = 350;
MaxGapDur = 150;

SyllEdges = 0:BinSize:MaxSyllDur;
GapEdges = 0:BinSize:MaxGapDur;

SyllGapDurHistFigure = figure;

Legend = [];
% First for pre days
for i = 1:Parameters.NoPreDays,
    Legend{end+1} = ['Pre #', num2str(i)];
    % First for directed song syllable durations
    subplot(2, 2, 1);
    hold on;
    NonZeroDurs = find(cellfun(@length, Parameters.PreDirSyllDurs{i}));
    plot(SyllEdges, histc(cell2mat(Parameters.PreDirSyllDurs{i}(NonZeroDurs)), SyllEdges) * 100 / sum(histc(cell2mat(Parameters.PreDirSyllDurs{i}(NonZeroDurs)), SyllEdges)), [Colours(i), Symbols(1), '-']);
    
    % Next for directed song gap durations
    subplot(2, 2, 2);
    hold on;
    NonZeroDurs = find(cellfun(@length, Parameters.PreDirGapDurs{i}));
    plot(GapEdges, histc(cell2mat(Parameters.PreDirGapDurs{i}(NonZeroDurs)), GapEdges) * 100 / sum(histc(cell2mat(Parameters.PreDirGapDurs{i}(NonZeroDurs)), GapEdges)), [Colours(i), Symbols(1), '-']);

    % Next for undirected song syllable durations
    subplot(2, 2, 3);
    hold on;
    NonZeroDurs = find(cellfun(@length, Parameters.PreUnDirSyllDurs{i}));
    plot(SyllEdges, histc(cell2mat(Parameters.PreUnDirSyllDurs{i}(NonZeroDurs)), SyllEdges) * 100 / sum(histc(cell2mat(Parameters.PreUnDirSyllDurs{i}(NonZeroDurs)), SyllEdges)), [Colours(i), Symbols(1), '-']);
    
    % Next for undirected song gap durations
    subplot(2, 2, 4);
    hold on;
    NonZeroDurs = find(cellfun(@length, Parameters.PreUnDirGapDurs{i}));
    plot(GapEdges, histc(cell2mat(Parameters.PreUnDirGapDurs{i}(NonZeroDurs)), GapEdges) * 100 / sum(histc(cell2mat(Parameters.PreUnDirGapDurs{i}(NonZeroDurs)), GapEdges)), [Colours(i), Symbols(1), '-']);
end

% Next for post days
for i = 1:Parameters.NoPostDays,
    Legend{end+1} = ['Post #', num2str(i)];
    % First for directed song syllable durations
    subplot(2, 2, 1);
    hold on;
    NonZeroDurs = find(cellfun(@length, Parameters.PostDirSyllDurs{i}));
    plot(SyllEdges, histc(cell2mat(Parameters.PostDirSyllDurs{i}(NonZeroDurs)), SyllEdges) * 100 / sum(histc(cell2mat(Parameters.PostDirSyllDurs{i}(NonZeroDurs)), SyllEdges)), [Colours(i), Symbols(2), '-']);
    
    % Next for directed song gap durations
    subplot(2, 2, 2);
    hold on;
    NonZeroDurs = find(cellfun(@length, Parameters.PostDirGapDurs{i}));
    plot(GapEdges, histc(cell2mat(Parameters.PostDirGapDurs{i}(NonZeroDurs)), GapEdges) * 100 / sum(histc(cell2mat(Parameters.PostDirGapDurs{i}(NonZeroDurs)), GapEdges)), [Colours(i), Symbols(2), '-']);

    % Next for undirected song syllable durations
    subplot(2, 2, 3);
    hold on;
    NonZeroDurs = find(cellfun(@length, Parameters.PostUnDirSyllDurs{i}));
    plot(SyllEdges, histc(cell2mat(Parameters.PostUnDirSyllDurs{i}(NonZeroDurs)), SyllEdges) * 100 / sum(histc(cell2mat(Parameters.PostUnDirSyllDurs{i}(NonZeroDurs)), SyllEdges)), [Colours(i), Symbols(2), '-']);
    
    % Next for undirected song gap durations
    subplot(2, 2, 4);
    hold on;
    NonZeroDurs = find(cellfun(@length, Parameters.PostUnDirGapDurs{i}));
    plot(GapEdges, histc(cell2mat(Parameters.PostUnDirGapDurs{i}(NonZeroDurs)), GapEdges) * 100 / sum(histc(cell2mat(Parameters.PostUnDirGapDurs{i}(NonZeroDurs)), GapEdges)), [Colours(i), Symbols(2), '-']);
end

figure(SyllGapDurHistFigure);
set(gcf, 'Color', 'k', 'Position', [400 150 1000 600]);
annotation('textbox', [0.45 0.95 0.1 0.05], 'String', Parameters.BirdName, 'Color', 'w', 'FontSize', 20)

for i = 1:4,
    subplot(2, 2, i);
    set(gca, 'Color', 'k');
    set(gca, 'XColor', 'w');
    set(gca, 'YColor', 'w');
    axis tight;
    Axis(i,:) = axis;
end
Axis(:,[1 3]) = 0;
Axis([1 3], 2) = MaxSyllDur;
Axis([2 4], 2) = MaxGapDur;
Axis([1 3], 4) = max(Axis([1 3], 4));
Axis([2 4], 4) = max(Axis([2 4], 4));

Title{1} = 'Dir - syllable durations';
Title{2} = 'Dir - gap durations';
Title{3} = 'Undir - syllable durations';
Title{4} = 'Undir - gap durations';

for i = 1:4,
    subplot(2, 2, i);
    legend(Legend, 'TextColor', 'w', 'FontSize', 14);
    set(gca, 'Box', 'on');
    axis(Axis(i,:));
    if (i > 2)
        xlabel('Duration (ms)', 'FontSize', 14);
    end
    
    if (mod(i, 2) == 1)
        ylabel('%', 'FontSize', 14);
    end
    title(Title{i}, 'FontSize', 14, 'Color', 'w');
end
%==========================================================================

%============= Identify bouts =============================================
InterBoutInterval = 500;

disp('Identifying bouts ...');

% First for pre song
for i = 1:Parameters.NoPreDays,
    % First for dir song
    InterBoutIntervalCellArray = num2cell(ones(size(Parameters.PreDirGapDurs{i}))*InterBoutInterval);
    if (size(Parameters.PreDirGapDurs{i},1) > size(Parameters.PreDirGapDurs{i},2))
        IndexCellArray = num2cell([1:1:length(Parameters.PreDirGapDurs{i})]');
    else
        IndexCellArray = num2cell([1:1:length(Parameters.PreDirGapDurs{i})]);
    end
    [Parameters.PreDirBouts{i}] = cellfun(@MA_IdentifyBouts, Parameters.PreDirGapDurs{i}, Parameters.PreDirLens{i}, Parameters.PreDirFs{i}, Parameters.PreDirOnsets{i}, Parameters.PreDirOffsets{i}, InterBoutIntervalCellArray, IndexCellArray, 'UniformOutput', 0);  

    % Then for undir song
    InterBoutIntervalCellArray = num2cell(ones(size(Parameters.PreUnDirGapDurs{i}))*InterBoutInterval);
    if (size(Parameters.PreUnDirGapDurs{i},1) > size(Parameters.PreUnDirGapDurs{i},2))
        IndexCellArray = num2cell([1:1:length(Parameters.PreUnDirGapDurs{i})]');
    else
        IndexCellArray = num2cell([1:1:length(Parameters.PreUnDirGapDurs{i})]);
    end
    [Parameters.PreUnDirBouts{i}] = cellfun(@MA_IdentifyBouts, Parameters.PreUnDirGapDurs{i}, Parameters.PreUnDirLens{i}, Parameters.PreUnDirFs{i}, Parameters.PreUnDirOnsets{i}, Parameters.PreUnDirOffsets{i}, InterBoutIntervalCellArray, IndexCellArray, 'UniformOutput', 0);  
end

% Next for post song
for i = 1:Parameters.NoPostDays,
    % First for dir song
    InterBoutIntervalCellArray = num2cell(ones(size(Parameters.PostDirGapDurs{i}))*InterBoutInterval);
    if (size(Parameters.PostDirGapDurs{i},1) > size(Parameters.PostDirGapDurs{i},2))
        IndexCellArray = num2cell([1:1:length(Parameters.PostDirGapDurs{i})]');
    else
        IndexCellArray = num2cell([1:1:length(Parameters.PostDirGapDurs{i})]);
    end
    [Parameters.PostDirBouts{i}] = cellfun(@MA_IdentifyBouts, Parameters.PostDirGapDurs{i}, Parameters.PostDirLens{i}, Parameters.PostDirFs{i}, Parameters.PostDirOnsets{i}, Parameters.PostDirOffsets{i}, InterBoutIntervalCellArray, IndexCellArray, 'UniformOutput', 0);  

    % Then for undir song
    InterBoutIntervalCellArray = num2cell(ones(size(Parameters.PostUnDirGapDurs{i}))*InterBoutInterval);
    if (size(Parameters.PostUnDirGapDurs{i},1) > size(Parameters.PostUnDirGapDurs{i},2))
        IndexCellArray = num2cell([1:1:length(Parameters.PostUnDirGapDurs{i})]');
    else
        IndexCellArray = num2cell([1:1:length(Parameters.PostUnDirGapDurs{i})]);
    end
    [Parameters.PostUnDirBouts{i}] = cellfun(@MA_IdentifyBouts, Parameters.PostUnDirGapDurs{i}, Parameters.PostUnDirLens{i}, Parameters.PostUnDirFs{i}, Parameters.PostUnDirOnsets{i}, Parameters.PostUnDirOffsets{i}, InterBoutIntervalCellArray, IndexCellArray, 'UniformOutput', 0);  
end

%==========================================================================

% %============= Now do the spectral correlations ===========================
% % First pick 10 random bouts for each day and each condition with bout
% % length > 3 seconds.
% 
% WindowSize = 10; % window size in ms for calculating PSD
% StepSize = 1; % step in ms for calculating PSD
% BandWidth = 1.5; % bandwith for generating multi-taper spectrogram
% BoutSize = 3; % minimum length of bout to be considered in sec
% 
% disp('Doing spectral correlations ...');
% % First for pre song
% 
% % First for undir song
% 
% FileTypeCellArray = cellstr(char(ones(length(Parameters.PreUnDirBouts), 1) * double(Parameters.FileType)))';
% WindowSizeCellArray = num2cell(ones(1, length(Parameters.PreUnDirBouts)) * WindowSize);
% StepSizeCellArray = num2cell(ones(1, length(Parameters.PreUnDirBouts)) * StepSize);
% BandWidthCellArray = num2cell(ones(1, length(Parameters.PreUnDirBouts)) * BandWidth);
% BoutSizeCellArray = num2cell(ones(1, length(Parameters.PreUnDirBouts)) * BoutSize);
% 
% [Parameters.PreUnDirSpectralCorrs] = cellfun(@MA_DoSpectralCorrelations, Parameters.PreUnDirBouts, Parameters.PreUnDirSongFileNames, Parameters.PreDataDir, FileTypeCellArray, WindowSizeCellArray, StepSizeCellArray, BandWidthCellArray, BoutSizeCellArray, 'UniformOutput', 0);
% 
% % Now for post song
% % First for undir song
% 
% FileTypeCellArray = cellstr(char(ones(length(Parameters.PostUnDirBouts), 1) * double(Parameters.FileType)))';
% WindowSizeCellArray = num2cell(ones(1, length(Parameters.PostUnDirBouts)) * WindowSize);
% StepSizeCellArray = num2cell(ones(1, length(Parameters.PostUnDirBouts)) * StepSize);
% BandWidthCellArray = num2cell(ones(1, length(Parameters.PostUnDirBouts)) * BandWidth);
% BoutSizeCellArray = num2cell(ones(1, length(Parameters.PostUnDirBouts)) * BoutSize);
% 
% [Parameters.PostUnDirSpectralCorrs] = cellfun(@MA_DoSpectralCorrelations, Parameters.PostUnDirBouts, Parameters.PostUnDirSongFileNames, Parameters.PostDataDir, FileTypeCellArray, WindowSizeCellArray, StepSizeCellArray, BandWidthCellArray, BoutSizeCellArray, 'UniformOutput', 0);
% 
% % First for dir song
% 
% FileTypeCellArray = cellstr(char(ones(length(Parameters.PreDirBouts), 1) * double(Parameters.FileType)))';
% WindowSizeCellArray = num2cell(ones(1, length(Parameters.PreDirBouts)) * WindowSize);
% StepSizeCellArray = num2cell(ones(1, length(Parameters.PreDirBouts)) * StepSize);
% BandWidthCellArray = num2cell(ones(1, length(Parameters.PreDirBouts)) * BandWidth);
% BoutSizeCellArray = num2cell(ones(1, length(Parameters.PreDirBouts)) * BoutSize);
% 
% [Parameters.PreDirSpectralCorrs] = cellfun(@MA_DoSpectralCorrelations, Parameters.PreDirBouts, Parameters.PreDirSongFileNames, Parameters.PreDataDir, FileTypeCellArray, WindowSizeCellArray, StepSizeCellArray, BandWidthCellArray, BoutSizeCellArray, 'UniformOutput', 0);
% 
% % Now for post song
% % First for undir song
% 
% FileTypeCellArray = cellstr(char(ones(length(Parameters.PostDirBouts), 1) * double(Parameters.FileType)))';
% WindowSizeCellArray = num2cell(ones(1, length(Parameters.PostDirBouts)) * WindowSize);
% StepSizeCellArray = num2cell(ones(1, length(Parameters.PostDirBouts)) * StepSize);
% BandWidthCellArray = num2cell(ones(1, length(Parameters.PostDirBouts)) * BandWidth);
% BoutSizeCellArray = num2cell(ones(1, length(Parameters.PostDirBouts)) * BoutSize);
% 
% [Parameters.PostDirSpectralCorrs] = cellfun(@MA_DoSpectralCorrelations, Parameters.PostDirBouts, Parameters.PostDirSongFileNames, Parameters.PostDataDir, FileTypeCellArray, WindowSizeCellArray, StepSizeCellArray, BandWidthCellArray, BoutSizeCellArray, 'UniformOutput', 0);
% %==========================================================================

%============= Template Matching ==========================================
TemplateMatchOutputDir = '/home/raghav/MicrolesionAnalysisResults_MT_Templates/';
[MotifTemplateDir, MotifTemplateFileName, MotifTemplateExt] = fileparts(Parameters.MotifTemplateFileName);

TemplateMatchOutputDir = [TemplateMatchOutputDir, MotifTemplateFileName, MotifTemplateExt, '.TemplateMatchResults', FileSep];
if (~exist(TemplateMatchOutputDir, 'dir'))
    mkdir(TemplateMatchOutputDir);
end

disp('Loading motif template ...');
Parameters.MotifTemplate = load(Parameters.MotifTemplateFileName);

disp('Doing template matching ...');

% First for pre song
for i = 1:Parameters.NoPreDays,
    % First for dir song
    disp(['   Pre Day #', num2str(i), ' - directed song ...']); 
    FileTypeCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double(Parameters.FileType)));
    RawDataDirCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double(Parameters.PreDataDir{i})));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double(TemplateMatchOutputDir)));
    LabelCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double('Motif')));
    TemplateTypeCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double('Spectrogram')));
    MotifTemplateCellArray = cell(size(Parameters.PreDirSongFileNames{i}));
    for j = 1:length(MotifTemplateCellArray),
        MotifTemplateCellArray{j} = Parameters.MotifTemplate;
    end
    % using cellfun so that i iterate over each element of the cell array.
    % To use cellfun, all of the other inputs also have to be in the form
    % of cell arrays of the same length - so the previous three lines
    % convert file type, data dir and output dir - common parameters for
    % all of the files into cell arrays
    cellfun(@MA_TemplateMatch, RawDataDirCellArray, Parameters.PreDirSongFileNames{i}, FileTypeCellArray, MotifTemplateCellArray, LabelCellArray, OutputDirCellArray, TemplateTypeCellArray, 'UniformOutput', 0);
    
    % next undirected songs
    disp(['   Pre Day #', num2str(i), ' - undirected song ...']); 

    FileTypeCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double(Parameters.FileType)));
    RawDataDirCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double(Parameters.PreDataDir{i})));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double(TemplateMatchOutputDir)));
    LabelCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double('Motif')));
    TemplateTypeCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double('Spectrogram')));
    MotifTemplateCellArray = cell(size(Parameters.PreUnDirSongFileNames{i}));
    for j = 1:length(MotifTemplateCellArray),
        MotifTemplateCellArray{j} = Parameters.MotifTemplate;
    end
    % using cellfun so that i iterate over each element of the cell array.
    % To use cellfun, all of the other inputs also have to be in the form
    % of cell arrays of the same length - so the previous three lines
    % convert file type, data dir and output dir - common parameters for
    % all of the files into cell arrays
    cellfun(@MA_TemplateMatch, RawDataDirCellArray, Parameters.PreUnDirSongFileNames{i}, FileTypeCellArray, MotifTemplateCellArray, LabelCellArray, OutputDirCellArray, TemplateTypeCellArray, 'UniformOutput', 0);
end

% Next for post song
for i = 1:Parameters.NoPostDays,
    % First for dir song
    disp(['   Post Day #', num2str(i), ' - directed song ...']); 
    FileTypeCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double(Parameters.FileType)));
    RawDataDirCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double(Parameters.PostDataDir{i})));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double(TemplateMatchOutputDir)));
    LabelCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double('Motif')));
    TemplateTypeCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double('Spectrogram')));
    MotifTemplateCellArray = cell(size(Parameters.PostDirSongFileNames{i}));
    for j = 1:length(MotifTemplateCellArray),
        MotifTemplateCellArray{j} = Parameters.MotifTemplate;
    end
    % using cellfun so that i iterate over each element of the cell array.
    % To use cellfun, all of the other inputs also have to be in the form
    % of cell arrays of the same length - so the previous three lines
    % convert file type, data dir and output dir - common parameters for
    % all of the files into cell arrays
    cellfun(@MA_TemplateMatch, RawDataDirCellArray, Parameters.PostDirSongFileNames{i}, FileTypeCellArray, MotifTemplateCellArray, LabelCellArray, OutputDirCellArray, TemplateTypeCellArray,  'UniformOutput', 0);
    
    % next undirected songs
    disp(['   Post Day #', num2str(i), ' - undirected song ...']); 

    FileTypeCellArray = cellstr(char(ones(length(Parameters.PostUnDirSongFileNames{i}), 1)*double(Parameters.FileType)));
    RawDataDirCellArray = cellstr(char(ones(length(Parameters.PostUnDirSongFileNames{i}), 1)*double(Parameters.PostDataDir{i})));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PostUnDirSongFileNames{i}), 1)*double(TemplateMatchOutputDir)));
    LabelCellArray = cellstr(char(ones(length(Parameters.PostUnDirSongFileNames{i}), 1)*double('Motif')));
    TemplateTypeCellArray = cellstr(char(ones(length(Parameters.PostUnDirSongFileNames{i}), 1)*double('Spectrogram')));
    MotifTemplateCellArray = cell(size(Parameters.PostUnDirSongFileNames{i}));
    for j = 1:length(MotifTemplateCellArray),
        MotifTemplateCellArray{j} = Parameters.MotifTemplate;
    end
    % using cellfun so that i iterate over each element of the cell array.
    % To use cellfun, all of the other inputs also have to be in the form
    % of cell arrays of the same length - so the previous three lines
    % convert file type, data dir and output dir - common parameters for
    % all of the files into cell arrays
    cellfun(@MA_TemplateMatch, RawDataDirCellArray, Parameters.PostUnDirSongFileNames{i}, FileTypeCellArray, MotifTemplateCellArray, LabelCellArray, OutputDirCellArray, TemplateTypeCellArray,  'UniformOutput', 0);
end
%==========================================================================

%============= Template Matching with other bird's songs ==================
TemplateMatchOutputDir = '/home/raghav/MicrolesionAnalysisResults_MT_Templates/OtherBirdSongs/';
[MotifTemplateDir, MotifTemplateFileName, MotifTemplateExt] = fileparts(Parameters.MotifTemplateFileName);

TemplateMatchOutputDir = [TemplateMatchOutputDir, MotifTemplateFileName, MotifTemplateExt, '.TemplateMatchResults', FileSep];
if (~exist(TemplateMatchOutputDir, 'dir'))
    mkdir(TemplateMatchOutputDir);
end

disp('Doing template matching to other birds songs...');

% First for wav files
for i = 1:Parameters.OtherBirdsSongDetails.NoPreDays,
    % undirected songs
    disp(['   Other birds songs ...']); 
    if (~isempty(Parameters.OtherBirdsSongDetails.PreUnDirSongFileNames{i}))
        FileTypeCellArray = cellstr(char(ones(length(Parameters.OtherBirdsSongDetails.PreUnDirSongFileNames{i}), 1)*double(Parameters.OtherBirdsSongDetails.PreDate{i})));
        RawDataDirCellArray = cellstr(char(ones(length(Parameters.OtherBirdsSongDetails.PreUnDirSongFileNames{i}), 1)*double(Parameters.OtherBirdsSongDetails.PreDataDir{i})));
        OutputDirCellArray = cellstr(char(ones(length(Parameters.OtherBirdsSongDetails.PreUnDirSongFileNames{i}), 1)*double(TemplateMatchOutputDir)));
        LabelCellArray = cellstr(char(ones(length(Parameters.OtherBirdsSongDetails.PreUnDirSongFileNames{i}), 1)*double('Motif')));
        TemplateTypeCellArray = cellstr(char(ones(length(Parameters.OtherBirdsSongDetails.PreUnDirSongFileNames{i}), 1)*double('Spectrogram')));
        MotifTemplateCellArray = cell(size(Parameters.OtherBirdsSongDetails.PreUnDirSongFileNames{i}));
        for j = 1:length(MotifTemplateCellArray),
            MotifTemplateCellArray{j} = Parameters.MotifTemplate;
        end
        % using cellfun so that i iterate over each element of the cell array.
        % To use cellfun, all of the other inputs also have to be in the form
        % of cell arrays of the same length - so the previous three lines
        % convert file type, data dir and output dir - common parameters for
        % all of the files into cell arrays
        cellfun(@MA_TemplateMatch, RawDataDirCellArray, Parameters.OtherBirdsSongDetails.PreUnDirSongFileNames{i}, FileTypeCellArray, MotifTemplateCellArray, LabelCellArray, OutputDirCellArray, TemplateTypeCellArray, 'UniformOutput', 0);
    end
end
%==========================================================================

%============= Shuffled song comparisons with template ====================
TemplateMatchOutputDir = '/home/raghav/MicrolesionAnalysisResults_MT_Templates/ShuffledSongComparisons/';
[MotifTemplateDir, MotifTemplateFileName, MotifTemplateExt] = fileparts(Parameters.MotifTemplateFileName);

TemplateMatchOutputDir = [TemplateMatchOutputDir, MotifTemplateFileName, MotifTemplateExt, '.TemplateMatchResults', FileSep];
if (~exist(TemplateMatchOutputDir, 'dir'))
    mkdir(TemplateMatchOutputDir);
end

disp('Doing shuffled song template matching ...');

% Only for one pre undir song day

disp(['   Pre Day #', num2str(i), ' - undirected song ...']); 

NumberofFiles = min(50, length(Parameters.PreUnDirSongFileNames{1}));
RandomFiles = Parameters.PreUnDirSongFileNames{1}(1:NumberofFiles);

FileTypeCellArray = cellstr(char(ones(length(RandomFiles), 1)*double(Parameters.FileType)));
RawDataDirCellArray = cellstr(char(ones(length(RandomFiles), 1)*double(Parameters.PreDataDir{1})));
OutputDirCellArray = cellstr(char(ones(length(RandomFiles), 1)*double(TemplateMatchOutputDir)));
LabelCellArray = cellstr(char(ones(length(RandomFiles), 1)*double('Motif')));
TemplateTypeCellArray = cellstr(char(ones(length(RandomFiles), 1)*double('Spectrogram')));
MotifTemplateCellArray = cell(size(RandomFiles));
for j = 1:length(MotifTemplateCellArray),
    MotifTemplateCellArray{j} = Parameters.MotifTemplate;
end
% using cellfun so that i iterate over each element of the cell array.
% To use cellfun, all of the other inputs also have to be in the form
% of cell arrays of the same length - so the previous three lines
% convert file type, data dir and output dir - common parameters for
% all of the files into cell arrays
cellfun(@MA_RandomTemplateMatch, RawDataDirCellArray, RandomFiles, FileTypeCellArray, MotifTemplateCellArray, LabelCellArray, OutputDirCellArray, TemplateTypeCellArray, 'UniformOutput', 0);
%==========================================================================

%============= Shuffled song comparisons with template shuffling only part of the song ====================
ShufflePercentages = [25 50 75];

for Shuffle = ShufflePercentages,
    TemplateMatchOutputDir = '/home/raghav/MicrolesionAnalysisResults_MT_Templates/PartShuffledSongComparisons/';
    [MotifTemplateDir, MotifTemplateFileName, MotifTemplateExt] = fileparts(Parameters.MotifTemplateFileName);

    TemplateMatchOutputDir = [TemplateMatchOutputDir, MotifTemplateFileName, MotifTemplateExt, '.Shuffle.', num2str(Shuffle), 'percent.TemplateMatchResults', FileSep];
    if (~exist(TemplateMatchOutputDir, 'dir'))
        mkdir(TemplateMatchOutputDir);
    end

    disp(['Doing ', num2str(Shuffle), ' % shuffled song template matching ...']);

    % Only for one pre undir song day

    disp(['   Pre Day #', num2str(i), ' - undirected song ...']); 

    NumberofFiles = min(50, length(Parameters.PreUnDirSongFileNames{1}));
    RandomFiles = Parameters.PreUnDirSongFileNames{1}(1:NumberofFiles);

    ShuffleCellArray = num2cell(ones(length(RandomFiles), 1) * Shuffle/100);
    FileTypeCellArray = cellstr(char(ones(length(RandomFiles), 1)*double(Parameters.FileType)));
    RawDataDirCellArray = cellstr(char(ones(length(RandomFiles), 1)*double(Parameters.PreDataDir{1})));
    OutputDirCellArray = cellstr(char(ones(length(RandomFiles), 1)*double(TemplateMatchOutputDir)));
    LabelCellArray = cellstr(char(ones(length(RandomFiles), 1)*double('Motif')));
    TemplateTypeCellArray = cellstr(char(ones(length(RandomFiles), 1)*double('Spectrogram')));
    MotifTemplateCellArray = cell(size(RandomFiles));
    for j = 1:length(MotifTemplateCellArray),
        MotifTemplateCellArray{j} = Parameters.MotifTemplate;
    end
    % using cellfun so that i iterate over each element of the cell array.
    % To use cellfun, all of the other inputs also have to be in the form
    % of cell arrays of the same length - so the previous three lines
    % convert file type, data dir and output dir - common parameters for
    % all of the files into cell arrays
    cellfun(@MA_PartRandomTemplateMatch, RawDataDirCellArray, RandomFiles, FileTypeCellArray, MotifTemplateCellArray, LabelCellArray, OutputDirCellArray, TemplateTypeCellArray, ShuffleCellArray, 'UniformOutput', 0);
end
%==========================================================================

%============= Shuffled song comparisons 1-100% for only one file of song ====================
% ShufflePercentages = [1:1:100];
% 
% for Shuffle = ShufflePercentages,
%     TemplateMatchOutputDir = '/home/raghav/MicrolesionAnalysisResults_MT_Templates/PartShuffledSongComparisons_1_to_100/';
%     [MotifTemplateDir, MotifTemplateFileName, MotifTemplateExt] = fileparts(Parameters.MotifTemplateFileName);
% 
%     TemplateMatchOutputDir = [TemplateMatchOutputDir, MotifTemplateFileName, MotifTemplateExt, '.Shuffle.', num2str(Shuffle), 'percent.TemplateMatchResults', FileSep];
%     if (~exist(TemplateMatchOutputDir, 'dir'))
%         mkdir(TemplateMatchOutputDir);
%     end
% 
%     disp(['Doing ', num2str(Shuffle), ' % shuffled song template matching ...']);
% 
%     % Only for one pre undir song day
% 
%     disp(['   Pre Day #', num2str(i), ' - undirected song ...']); 
% 
%     NumberofFiles = min(2, length(Parameters.PreUnDirSongFileNames{1}));
%     RandomFiles = Parameters.PreUnDirSongFileNames{1}(1:NumberofFiles);
% 
%     ShuffleCellArray = num2cell(ones(length(RandomFiles), 1) * Shuffle/100);
%     FileTypeCellArray = cellstr(char(ones(length(RandomFiles), 1)*double(Parameters.FileType)));
%     RawDataDirCellArray = cellstr(char(ones(length(RandomFiles), 1)*double(Parameters.PreDataDir{1})));
%     OutputDirCellArray = cellstr(char(ones(length(RandomFiles), 1)*double(TemplateMatchOutputDir)));
%     LabelCellArray = cellstr(char(ones(length(RandomFiles), 1)*double('Motif')));
%     TemplateTypeCellArray = cellstr(char(ones(length(RandomFiles), 1)*double('Spectrogram')));
%     MotifTemplateCellArray = cell(size(RandomFiles));
%     for j = 1:length(MotifTemplateCellArray),
%         MotifTemplateCellArray{j} = Parameters.MotifTemplate;
%     end
%     % using cellfun so that i iterate over each element of the cell array.
%     % To use cellfun, all of the other inputs also have to be in the form
%     % of cell arrays of the same length - so the previous three lines
%     % convert file type, data dir and output dir - common parameters for
%     % all of the files into cell arrays
%     cellfun(@MA_PartRandomTemplateMatch, RawDataDirCellArray, RandomFiles, FileTypeCellArray, MotifTemplateCellArray, LabelCellArray, OutputDirCellArray, TemplateTypeCellArray, ShuffleCellArray, 'UniformOutput', 0);
% end
%==========================================================================

% %============= Shuffled song comparisons (time bins only) 1-100% for only one file of song ====================
% ShufflePercentages = 100;
% 
% for Shuffle = ShufflePercentages,
%     TemplateMatchOutputDir = '/home/raghav/MicrolesionAnalysisResults_MT_Templates/PartTimeShuffledSongComparisons_1_to_100/';
%     [MotifTemplateDir, MotifTemplateFileName, MotifTemplateExt] = fileparts(Parameters.MotifTemplateFileName);
% 
%     TemplateMatchOutputDir = [TemplateMatchOutputDir, MotifTemplateFileName, MotifTemplateExt, '.Shuffle.', num2str(Shuffle), 'percent.TemplateMatchResults', FileSep];
%     if (~exist(TemplateMatchOutputDir, 'dir'))
%         mkdir(TemplateMatchOutputDir);
%     end
% 
%     disp(['Doing ', num2str(Shuffle), ' % shuffled song template matching ...']);
% 
%     % Only for one pre undir song day
% 
%     disp(['   Pre Day #', num2str(i), ' - undirected song ...']); 
% 
%     NumberofFiles = min(2, length(Parameters.PreUnDirSongFileNames{1}));
%     RandomFiles = Parameters.PreUnDirSongFileNames{1}(1:NumberofFiles);
% 
%     ShuffleCellArray = num2cell(ones(length(RandomFiles), 1) * Shuffle/100);
%     FileTypeCellArray = cellstr(char(ones(length(RandomFiles), 1)*double(Parameters.FileType)));
%     RawDataDirCellArray = cellstr(char(ones(length(RandomFiles), 1)*double(Parameters.PreDataDir{1})));
%     OutputDirCellArray = cellstr(char(ones(length(RandomFiles), 1)*double(TemplateMatchOutputDir)));
%     LabelCellArray = cellstr(char(ones(length(RandomFiles), 1)*double('Motif')));
%     TemplateTypeCellArray = cellstr(char(ones(length(RandomFiles), 1)*double('Spectrogram')));
%     MotifTemplateCellArray = cell(size(RandomFiles));
%     for j = 1:length(MotifTemplateCellArray),
%         MotifTemplateCellArray{j} = Parameters.MotifTemplate;
%     end
%     % using cellfun so that i iterate over each element of the cell array.
%     % To use cellfun, all of the other inputs also have to be in the form
%     % of cell arrays of the same length - so the previous three lines
%     % convert file type, data dir and output dir - common parameters for
%     % all of the files into cell arrays
%     cellfun(@MA_PartRandomTimeTemplateMatch, RawDataDirCellArray, RandomFiles, FileTypeCellArray, MotifTemplateCellArray, LabelCellArray, OutputDirCellArray, TemplateTypeCellArray, ShuffleCellArray, 'UniformOutput', 0);
% end
% %==========================================================================

%============= Load results of template matching ==========================
% Next for the 100% shuffled song comparisons
TemplateMatchOutputDir = '/home/raghav/MicrolesionAnalysisResults_MT_Templates/ShuffledSongComparisons/';
[MotifTemplateDir, MotifTemplateFileName, MotifTemplateExt] = fileparts(Parameters.MotifTemplateFileName);

TemplateMatchOutputDir = [TemplateMatchOutputDir, MotifTemplateFileName, MotifTemplateExt, '.TemplateMatchResults', FileSep];
disp('Loading results of shuffled song template matching ...');

NumberofFiles = min(50, length(Parameters.PreUnDirSongFileNames{1}));
RandomFiles = Parameters.PreUnDirSongFileNames{1}(1:NumberofFiles);

LabelCellArray = cellstr(char(ones(length(RandomFiles), 1)*double('Motif')));
OutputDirCellArray = cellstr(char(ones(length(RandomFiles), 1)*double(TemplateMatchOutputDir)));

Parameters.ShuffledSongMatches = cellfun(@MA_LoadShuffledSong_TemplateMatchResultsFile, RandomFiles, OutputDirCellArray, LabelCellArray, 'UniformOutput', 0);
TempShuffledSongMatches = [];
TempMaxShuffledSongMatches = [];
for i = 1:length(Parameters.ShuffledSongMatches),
    TempMaxShuffledSongMatches = [TempMaxShuffledSongMatches; cell2mat(cellfun(@max, Parameters.ShuffledSongMatches{i}, 'UniformOutput', 0))'];
    TempShuffledSongMatches = [TempShuffledSongMatches; cell2mat(Parameters.ShuffledSongMatches{i}(:))];
end
Parameters.MaxShuffledSongMatches = TempMaxShuffledSongMatches;
Parameters.ShuffledSongMatches = TempShuffledSongMatches;
clear TempShuffledSongMatches TempMaxShuffledSongMatches;

% ShufflePercentages = [1:1:100];
% 
% ShuffleIndex = 0;
% for Shuffle = ShufflePercentages,
% 
%     ShuffleIndex = ShuffleIndex + 1;
%     
%     TemplateMatchOutputDir = '/home/raghav/MicrolesionAnalysisResults_MT_Templates/PartShuffledSongComparisons_1_to_100/';
%     [MotifTemplateDir, MotifTemplateFileName, MotifTemplateExt] = fileparts(Parameters.MotifTemplateFileName);
% 
%     TemplateMatchOutputDir = [TemplateMatchOutputDir, MotifTemplateFileName, MotifTemplateExt, '.Shuffle.', num2str(Shuffle), 'percent.TemplateMatchResults', FileSep];
% 
%     disp(['Loading results of ', num2str(Shuffle), ' % shuffled song template matching ...']);
% 
%     NumberofFiles = min(2, length(Parameters.PreUnDirSongFileNames{1}));
%     RandomFiles = Parameters.PreUnDirSongFileNames{1}(1:NumberofFiles);
% 
%     LabelCellArray = cellstr(char(ones(length(RandomFiles), 1)*double('Motif')));
%     OutputDirCellArray = cellstr(char(ones(length(RandomFiles), 1)*double(TemplateMatchOutputDir)));
% 
%     Parameters.PartShuffledSongMatchesAll{ShuffleIndex} = cellfun(@MA_LoadShuffledSong_TemplateMatchResultsFile, RandomFiles, OutputDirCellArray, LabelCellArray, 'UniformOutput', 0);
%     Parameters.PartShuffledSongMatchesAll{ShuffleIndex} = Parameters.PartShuffledSongMatchesAll{ShuffleIndex}(:)';
%     Parameters.PartShuffledSongMatchesAll{ShuffleIndex} = cell2mat(Parameters.PartShuffledSongMatchesAll{ShuffleIndex});
%     Parameters.PartShuffledSongMatchesAll{ShuffleIndex} = Parameters.PartShuffledSongMatchesAll{ShuffleIndex}(:);
%     Parameters.PartShuffleTemplateMatchThresholdAll{ShuffleIndex} = mean(Parameters.PartShuffledSongMatchesAll{ShuffleIndex}) + 3*std(Parameters.PartShuffledSongMatchesAll{ShuffleIndex});
% end
% 
% %ShuffledMatchesFigure = figure;
% figure(ShuffledMatchesFigure);
% hold on;
% TemplateMatchOutputDir = '/home/raghav/MicrolesionAnalysisResults_MT_Templates/';
% [MotifTemplateDir, MotifTemplateFileName, MotifTemplateExt] = fileparts(Parameters.MotifTemplateFileName);
% TemplateMatchOutputDir = [TemplateMatchOutputDir, MotifTemplateFileName, MotifTemplateExt, '.TemplateMatchResults', FileSep];
% TempNormalSongFile = load([TemplateMatchOutputDir, Parameters.PreUnDirSongFileNames{1}{1}, '.Motif.TempMatch.mat']);
% TempPartShuffledMatches = cell2mat(Parameters.PartShuffledSongMatchesAll);
% plot(mean(TempPartShuffledMatches(1:2,:)), 'ko-');

% ShufflePercentages = 100;
% 
% ShuffleIndex = 0;
% for Shuffle = ShufflePercentages,
% 
%     ShuffleIndex = ShuffleIndex + 1;
%     
%     TemplateMatchOutputDir = '/home/raghav/MicrolesionAnalysisResults_MT_Templates/PartTimeShuffledSongComparisons_1_to_100/';
%     [MotifTemplateDir, MotifTemplateFileName, MotifTemplateExt] = fileparts(Parameters.MotifTemplateFileName);
% 
%     TemplateMatchOutputDir = [TemplateMatchOutputDir, MotifTemplateFileName, MotifTemplateExt, '.Shuffle.', num2str(Shuffle), 'percent.TemplateMatchResults', FileSep];
% 
%     disp(['Loading results of ', num2str(Shuffle), ' % shuffled song template matching ...']);
% 
%     NumberofFiles = min(2, length(Parameters.PreUnDirSongFileNames{1}));
%     RandomFiles = Parameters.PreUnDirSongFileNames{1}(1:NumberofFiles);
% 
%     LabelCellArray = cellstr(char(ones(length(RandomFiles), 1)*double('Motif')));
%     OutputDirCellArray = cellstr(char(ones(length(RandomFiles), 1)*double(TemplateMatchOutputDir)));
% 
%     Parameters.PartTimeShuffledSongMatchesAll{ShuffleIndex} = cellfun(@MA_LoadShuffledSong_TemplateMatchResultsFile, RandomFiles, OutputDirCellArray, LabelCellArray, 'UniformOutput', 0);
%     Parameters.PartTimeShuffledSongMatchesAll{ShuffleIndex} = Parameters.PartTimeShuffledSongMatchesAll{ShuffleIndex}(:)';
%     Parameters.PartTimeShuffledSongMatchesAll{ShuffleIndex} = cell2mat(Parameters.PartTimeShuffledSongMatchesAll{ShuffleIndex});
%     Parameters.PartTimeShuffledSongMatchesAll{ShuffleIndex} = Parameters.PartTimeShuffledSongMatchesAll{ShuffleIndex}(:);
%     Parameters.PartTimeShuffleTemplateMatchThresholdAll{ShuffleIndex} = mean(Parameters.PartTimeShuffledSongMatchesAll{ShuffleIndex}) + 3*std(Parameters.PartTimeShuffledSongMatchesAll{ShuffleIndex});
% end
% 
% ShuffledMatchesFigure = figure;
% figure(ShuffledMatchesFigure);
% hold on;
% TemplateMatchOutputDir = '/home/raghav/MicrolesionAnalysisResults_MT_Templates/';
% [MotifTemplateDir, MotifTemplateFileName, MotifTemplateExt] = fileparts(Parameters.MotifTemplateFileName);
% TemplateMatchOutputDir = [TemplateMatchOutputDir, MotifTemplateFileName, MotifTemplateExt, '.TemplateMatchResults', FileSep];
% TempNormalSongFile = load([TemplateMatchOutputDir, Parameters.PreUnDirSongFileNames{1}{1}, '.Motif.TempMatch.mat']);
% TempPartTimeShuffledMatches = cell2mat(Parameters.PartTimeShuffledSongMatchesAll);
% plot(mean(TempPartTimeShuffledMatches(1:2,:)), 'ko-');

ShufflePercentages = [25 50 75];

ShuffleIndex = 0;
for Shuffle = ShufflePercentages,

    ShuffleIndex = ShuffleIndex + 1;
    
    TemplateMatchOutputDir = '/home/raghav/MicrolesionAnalysisResults_MT_Templates/PartShuffledSongComparisons/';
    [MotifTemplateDir, MotifTemplateFileName, MotifTemplateExt] = fileparts(Parameters.MotifTemplateFileName);

    TemplateMatchOutputDir = [TemplateMatchOutputDir, MotifTemplateFileName, MotifTemplateExt, '.Shuffle.', num2str(Shuffle), 'percent.TemplateMatchResults', FileSep];

    disp(['Loading results of ', num2str(Shuffle), ' % shuffled song template matching ...']);

    NumberofFiles = min(50, length(Parameters.PreUnDirSongFileNames{1}));
    RandomFiles = Parameters.PreUnDirSongFileNames{1}(1:NumberofFiles);

    LabelCellArray = cellstr(char(ones(length(RandomFiles), 1)*double('Motif')));
    OutputDirCellArray = cellstr(char(ones(length(RandomFiles), 1)*double(TemplateMatchOutputDir)));

    Parameters.PartShuffledSongMatches{ShuffleIndex} = cellfun(@MA_LoadShuffledSong_TemplateMatchResultsFile, RandomFiles, OutputDirCellArray, LabelCellArray, 'UniformOutput', 0);
    
    TempShuffledSongMatches = [];
    TempMaxShuffledSongMatches = [];
    for i = 1:length(Parameters.PartShuffledSongMatches{ShuffleIndex}),
        TempMaxShuffledSongMatches = [TempMaxShuffledSongMatches; cell2mat(cellfun(@max, Parameters.PartShuffledSongMatches{ShuffleIndex}{i}, 'UniformOutput', 0))'];
        TempShuffledSongMatches = [TempShuffledSongMatches; cell2mat(Parameters.PartShuffledSongMatches{ShuffleIndex}{i}(:))];
    end
    Parameters.PartMaxShuffledSongMatches{ShuffleIndex} = TempMaxShuffledSongMatches;
    Parameters.PartShuffledSongMatches{ShuffleIndex} = TempShuffledSongMatches;
    clear TempShuffledSongMatches TempMaxShuffledSongMatches;
end

ActualTemplateMatchThreshold = max(Parameters.PartMaxShuffledSongMatches{2}); % 50% shuffled songs is used as thresholdls 
Results.ActualTemplateMatchThreshold = ActualTemplateMatchThreshold;

% Loading results for other birds songs

TemplateMatchOutputDir = '/home/raghav/MicrolesionAnalysisResults_MT_Templates/OtherBirdSongs/';
[MotifTemplateDir, MotifTemplateFileName, MotifTemplateExt] = fileparts(Parameters.MotifTemplateFileName);

TemplateMatchOutputDir = [TemplateMatchOutputDir, MotifTemplateFileName, MotifTemplateExt, '.TemplateMatchResults', FileSep];

disp('Loading results of other bird song template matching ...');

% Load results of other birds songs
for i = 1:Parameters.OtherBirdsSongDetails.NoPreDays,
    % undirected songs
    ZeroThreshold = 0;
    disp(['   Pre Day #', num2str(i), ' - other bird song ...']); 
    if (~isempty(Parameters.OtherBirdsSongDetails.PreUnDirSongFileNames{i}))
        LabelCellArray = cellstr(char(ones(length(Parameters.OtherBirdsSongDetails.PreUnDirSongFileNames{i}), 1)*double('Motif')));
        OutputDirCellArray = cellstr(char(ones(length(Parameters.OtherBirdsSongDetails.PreUnDirSongFileNames{i}), 1)*double(TemplateMatchOutputDir)));
        ThresholdCellArray = num2cell(ones(length(Parameters.OtherBirdsSongDetails.PreUnDirSongFileNames{i}), 1) * ZeroThreshold);
        SongFileNoCellArray = num2cell([1:1:length(Parameters.OtherBirdsSongDetails.PreUnDirSongFileNames{i})]');
        % using cellfun so that i iterate over each element of the cell array.
        % To use cellfun, all of the other inputs also have to be in the form
        % of cell arrays of the same length - so the previous three lines
        % convert file type, data dir and output dir - common parameters for
        % all of the files into cell arrays
        Parameters.OtherBirdSongResults{i} = cellfun(@MA_LoadTemplateMatchResultsFile, Parameters.OtherBirdsSongDetails.PreUnDirSongFileNames{i}, OutputDirCellArray, ThresholdCellArray, LabelCellArray, SongFileNoCellArray, 'UniformOutput', 0);
    end
end

% Loading results for pre song with threshold as 0 - basically all the
% peaks
TemplateMatchOutputDir = '/home/raghav/MicrolesionAnalysisResults_MT_Templates/';
[MotifTemplateDir, MotifTemplateFileName, MotifTemplateExt] = fileparts(Parameters.MotifTemplateFileName);

TemplateMatchOutputDir = [TemplateMatchOutputDir, MotifTemplateFileName, MotifTemplateExt, '.TemplateMatchResults', FileSep];

disp('Loading results of song template matching ...');

% First for pre song - threshold 0, loading results for pre undir alone
for i = 1:Parameters.NoPreDays,
    % next undirected songs
    ZeroThreshold = 0;
    disp(['   Pre Day #', num2str(i), ' - undirected song ...']); 

    LabelCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double('Motif')));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double(TemplateMatchOutputDir)));
    ThresholdCellArray = num2cell(ones(length(Parameters.PreUnDirSongFileNames{i}), 1) * ZeroThreshold);
    SongFileNoCellArray = num2cell([1:1:length(Parameters.PreUnDirSongFileNames{i})]');
    % using cellfun so that i iterate over each element of the cell array.
    % To use cellfun, all of the other inputs also have to be in the form
    % of cell arrays of the same length - so the previous three lines
    % convert file type, data dir and output dir - common parameters for
    % all of the files into cell arrays
    Parameters.PreZeroThresholdUnDirResults{i} = cellfun(@MA_LoadTemplateMatchResultsFile, Parameters.PreUnDirSongFileNames{i}, OutputDirCellArray, ThresholdCellArray, LabelCellArray, SongFileNoCellArray, 'UniformOutput', 0);
end

disp('Loading results of song template matching ...');
% Load results for all days with zero threshold - can then modify threshold
% later
% First for pre song
for i = 1:Parameters.NoPreDays,
    % First for dir song
    disp(['   Pre Day #', num2str(i), ' - directed song ...']); 
    ZeroThreshold = 0;
    LabelCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double('Motif')));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double(TemplateMatchOutputDir)));
    ThresholdCellArray = num2cell(ones(length(Parameters.PreDirSongFileNames{i}), 1) * ZeroThreshold);
    SongFileNoCellArray = num2cell([1:1:length(Parameters.PreDirSongFileNames{i})]');
    % using cellfun so that i iterate over each element of the cell array.
    % To use cellfun, all of the other inputs also have to be in the form
    % of cell arrays of the same length - so the previous three lines
    % convert file type, data dir and output dir - common parameters for
    % all of the files into cell arrays
    Parameters.PreDirResults{i} = cellfun(@MA_LoadTemplateMatchResultsFile, Parameters.PreDirSongFileNames{i}, OutputDirCellArray, ThresholdCellArray, LabelCellArray, SongFileNoCellArray, 'UniformOutput', 0);
    
    % next undirected songs
    disp(['   Pre Day #', num2str(i), ' - undirected song ...']); 

    LabelCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double('Motif')));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double(TemplateMatchOutputDir)));
    ThresholdCellArray = num2cell(ones(length(Parameters.PreUnDirSongFileNames{i}), 1) * ZeroThreshold);
    SongFileNoCellArray = num2cell([1:1:length(Parameters.PreUnDirSongFileNames{i})]');
    % using cellfun so that i iterate over each element of the cell array.
    % To use cellfun, all of the other inputs also have to be in the form
    % of cell arrays of the same length - so the previous three lines
    % convert file type, data dir and output dir - common parameters for
    % all of the files into cell arrays
    Parameters.PreUnDirResults{i} = cellfun(@MA_LoadTemplateMatchResultsFile, Parameters.PreUnDirSongFileNames{i}, OutputDirCellArray, ThresholdCellArray, LabelCellArray, SongFileNoCellArray, 'UniformOutput', 0);
end

% Next for post song
for i = 1:Parameters.NoPostDays,
    ZeroThreshold = 0;
    % First for dir song
    disp(['   Post Day #', num2str(i), ' - directed song ...']); 
    LabelCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double('Motif')));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double(TemplateMatchOutputDir)));
    ThresholdCellArray = num2cell(ones(length(Parameters.PostDirSongFileNames{i}), 1) * ZeroThreshold);
    SongFileNoCellArray = num2cell([1:1:length(Parameters.PostDirSongFileNames{i})]');
    % using cellfun so that i iterate over each element of the cell array.
    % To use cellfun, all of the other inputs also have to be in the form
    % of cell arrays of the same length - so the previous three lines
    % convert file type, data dir and output dir - common parameters for
    % all of the files into cell arrays
    Parameters.PostDirResults{i} = cellfun(@MA_LoadTemplateMatchResultsFile, Parameters.PostDirSongFileNames{i}, OutputDirCellArray, ThresholdCellArray, LabelCellArray, SongFileNoCellArray, 'UniformOutput', 0);
    
    % next undirected songs
    disp(['   Post Day #', num2str(i), ' - undirected song ...']); 

    LabelCellArray = cellstr(char(ones(length(Parameters.PostUnDirSongFileNames{i}), 1)*double('Motif')));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PostUnDirSongFileNames{i}), 1)*double(TemplateMatchOutputDir)));
    ThresholdCellArray = num2cell(ones(length(Parameters.PostUnDirSongFileNames{i}), 1) * ZeroThreshold);
    SongFileNoCellArray = num2cell([1:1:length(Parameters.PostUnDirSongFileNames{i})]');
    % using cellfun so that i iterate over each element of the cell array.
    % To use cellfun, all of the other inputs also have to be in the form
    % of cell arrays of the same length - so the previous three lines
    % convert file type, data dir and output dir - common parameters for
    % all of the files into cell arrays
    Parameters.PostUnDirResults{i} = cellfun(@MA_LoadTemplateMatchResultsFile, Parameters.PostUnDirSongFileNames{i}, OutputDirCellArray, ThresholdCellArray, LabelCellArray, SongFileNoCellArray, 'UniformOutput', 0);
end
%==========================================================================


%============= Plot match results  ========================================
% disp('Plotting match results ...');
% figure;
% 
% Index = 1;
% for i = 1:length(Parameters.PreDirResults),
%     DirAllMatches{i} = cell2mat(Parameters.PreDirResults{i});
%     DirTotalSongTime{i} = cell2mat(Parameters.PreDirLens{i});
%     
%     if (~isempty(DirAllMatches{i}))
%         Results.DirAllMatches(Index,:) = [i mean(DirAllMatches{i}(:,1)) std(DirAllMatches{i}(:,1))];
%         Results.DirNumMatchesPerSec(Index,:) = [i size(DirAllMatches{i},1)*1000/sum(DirTotalSongTime{i})];
%         
% %         subplot(2,1,1);
% %         errorbar(i, mean(DirAllMatches{i}(:,1)), std(DirAllMatches{i}(:,1)), 'ro', 'MarkerSize', 6);
% %         hold on;
% % 
% %         subplot(2,1,2);
% %         plot(i, size(DirAllMatches{i},1)*1000/sum(DirTotalSongTime{i}), 'ro', 'MarkerSize', 6);
% %         hold on;
%     else
%         Results.DirAllMatches(Index,:) = [i NaN NaN];
%         Results.DirNumMatchesPerSec(Index,:) = [i 0];
%     end
%     
%     UnDirAllMatches{i} = cell2mat(Parameters.PreUnDirResults{i});
%     UnDirTotalSongTime{i} = cell2mat(Parameters.PreUnDirLens{i});
%     
%     if (~isempty(UnDirAllMatches{i}))
%         Results.UnDirAllMatches(Index,:) = [i mean(UnDirAllMatches{i}(:,1)) std(UnDirAllMatches{i}(:,1))];
%         Results.UnDirNumMatchesPerSec(Index,:) = [i size(UnDirAllMatches{i},1)*1000/sum(UnDirTotalSongTime{i})];
%         
% %         subplot(2,1,1);
% %         errorbar(i, mean(UnDirAllMatches{i}(:,1)), std(UnDirAllMatches{i}(:,1)), 'bo', 'MarkerSize', 6);
% %         hold on;
% % 
% %         subplot(2,1,2);
% %         plot(i, size(UnDirAllMatches{i},1)*1000/sum(UnDirTotalSongTime{i}), 'bo', 'MarkerSize', 6);
% %         hold on;
%     else
%         Results.UnDirAllMatches(Index,:) = [i NaN NaN];
%         Results.UnDirNumMatchesPerSec(Index,:) = [i 0];
%     end
%     Index = Index + 1;
% end
% 
% for i = 1:length(Parameters.PostDirResults),
%     DirAllMatches{i + length(Parameters.PreDirResults)} = cell2mat(Parameters.PostDirResults{i});
%     DirTotalSongTime{i + length(Parameters.PreDirResults)} = cell2mat(Parameters.PostDirLens{i});
%     
%     if (~isempty(DirAllMatches{i + length(Parameters.PreDirResults)}))
%         Results.DirAllMatches(Index,:) = [(i + length(Parameters.PreDirResults)) mean(DirAllMatches{i + length(Parameters.PreDirResults)}(:,1)) std(DirAllMatches{i + length(Parameters.PreDirResults)}(:,1))];
%         Results.DirNumMatchesPerSec(Index,:) = [(i + length(Parameters.PreDirResults)) size(DirAllMatches{i + length(Parameters.PreDirResults)},1)*1000/sum(DirTotalSongTime{i + length(Parameters.PreDirResults)})];
% 
% %         subplot(2,1,1);
% %         errorbar(i + length(Parameters.PreDirResults), mean(DirAllMatches{i + length(Parameters.PreDirResults)}(:,1)), std(DirAllMatches{i + length(Parameters.PreDirResults)}(:,1)), 'ro', 'MarkerSize', 6);
% %         hold on;
% % 
% %         subplot(2,1,2);
% %         plot(i + length(Parameters.PreDirResults), size(DirAllMatches{i + length(Parameters.PreDirResults)},1)*1000/sum(DirTotalSongTime{i + length(Parameters.PreDirResults)}), 'ro', 'MarkerSize', 6);
% %         hold on;
%     else
%         Results.DirAllMatches(Index,:) = [(i + length(Parameters.PreDirResults)) NaN NaN];
%         Results.DirNumMatchesPerSec(Index,:) = [(i + length(Parameters.PreDirResults)) 0];
%     end
%     
%     UnDirAllMatches{i + length(Parameters.PreDirResults)} = cell2mat(Parameters.PostUnDirResults{i});
%     UnDirTotalSongTime{i + length(Parameters.PreDirResults)} = cell2mat(Parameters.PostUnDirLens{i});
%     
%     if (~isempty(UnDirAllMatches{i + length(Parameters.PreDirResults)}))
%         Results.UnDirAllMatches(Index,:) = [(i + length(Parameters.PreDirResults)) mean(UnDirAllMatches{i + length(Parameters.PreDirResults)}(:,1)) std(UnDirAllMatches{i + length(Parameters.PreDirResults)}(:,1))];
%         Results.UnDirNumMatchesPerSec(Index,:) = [(i + length(Parameters.PreDirResults)) size(UnDirAllMatches{i + length(Parameters.PreDirResults)},1)*1000/sum(UnDirTotalSongTime{i + length(Parameters.PreDirResults)})];
% 
% %         subplot(2,1,1);
% %         errorbar(i + length(Parameters.PreDirResults), mean(UnDirAllMatches{i + length(Parameters.PreDirResults)}(:,1)), std(UnDirAllMatches{i + length(Parameters.PreDirResults)}(:,1)), 'bo', 'MarkerSize', 6);
% %         hold on;
% % 
% %         subplot(2,1,2);
% %         plot(i + length(Parameters.PreDirResults), size(UnDirAllMatches{i + length(Parameters.PreDirResults)},1)*1000/sum(UnDirTotalSongTime{i + length(Parameters.PreDirResults)}), 'bo', 'MarkerSize', 6);
% %         hold on;
%     else
%         Results.UnDirAllMatches(Index,:) = [(i + length(Parameters.PreDirResults)) NaN NaN];
%         Results.UnDirNumMatchesPerSec(Index,:) = [(i + length(Parameters.PreDirResults)) 0];
%     end
%     Index = Index + 1;
% end
% 
% subplot(2,1,1);
% errorbar(Results.DirAllMatches(:,1), Results.DirAllMatches(:,2), Results.DirAllMatches(:,3), 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
% hold on;
% errorbar(Results.UnDirAllMatches(:,1), Results.UnDirAllMatches(:,2), Results.UnDirAllMatches(:,3), 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
% 
% axis tight;
% Temp = axis;
% Temp = [0.5 (length(Parameters.PreDirResults) + length(Parameters.PostDirResults) + 0.5) 0 1.05*Temp(4)];
% axis(Temp);
% 
% plot([length(Parameters.PreDirResults) length(Parameters.PreDirResults)]+0.5, [0 Temp(4)], 'w--', 'LineWidth', 2);
% plot([Temp(1) Temp(2)], [ActualTemplateMatchThreshold ActualTemplateMatchThreshold], 'w--', 'LineWidth', 2);
% ylabel('Template match value', 'FontSize', 16);
% text((length(Parameters.PreDirResults) + 0.45), Temp(4)/20, ['Microlesion surgery: ', Parameters.SurgeryDate], 'Rotation', 90, 'Color', 'w', 'FontSize', 10);
% text(0.75, ActualTemplateMatchThreshold + Temp(4)/20, ['Template match threshold'], 'Color', 'w', 'FontSize', 10);
% 
% subplot(2,1,2);
% plot(Results.DirNumMatchesPerSec(:,1), Results.DirNumMatchesPerSec(:,2), 'ro-', 'MarkerSize', 8, 'LineWidth', 2);
% hold on;
% plot(Results.UnDirNumMatchesPerSec(:,1), Results.UnDirNumMatchesPerSec(:,2), 'bo-', 'MarkerSize', 8, 'LineWidth', 2);
% 
% axis tight;
% Temp = axis;
% Temp = [0.5 (length(Parameters.PreDirResults) + length(Parameters.PostDirResults) + 0.5) 0 1.05*Temp(4)];
% axis(Temp);
% plot([length(Parameters.PreDirResults) length(Parameters.PreDirResults)]+0.5, [0 Temp(4)], 'w--', 'LineWidth', 2);
% ylabel('No of matches / sec', 'FontSize', 16);
% text((length(Parameters.PreDirResults) + 0.45), Temp(4)/20, ['Microlesion surgery: ', Parameters.SurgeryDate], 'Rotation', 90, 'Color', 'w', 'FontSize', 10);
% 
% for i = 1:length(Parameters.PreDirResults),
%     XLabelString{i} = Parameters.PreDate{i};
% end
% 
% for i = 1:length(Parameters.PostDirResults),
%     XLabelString{i + length(Parameters.PreDirResults)} = Parameters.PostDate{i};
% end
% 
% set(gca, 'XTick', (1:1:(length(Parameters.PostDirResults) + length(Parameters.PreDirResults))), 'XTickLabel', XLabelString);
% set(gcf, 'Color', 'k', 'Position', [400 150 1000 600]);
% annotation('textbox', [0.45 0.95 0.1 0.05], 'String', Parameters.BirdName, 'Color', 'w', 'FontSize', 20)
% 
% subplot(2,1,1);
% set(gca, 'XTick', (1:1:(length(Parameters.PostDirResults) + length(Parameters.PreDirResults))), 'XTickLabel', XLabelString);
% 
% for i = 1:2,
%     subplot(2, 1, i);
%     set(gca, 'Color', 'k');
%     set(gca, 'XColor', 'w');
%     set(gca, 'YColor', 'w');
%     set(gca, 'FontSize', 14);
%     set(gca, 'Box', 'on');
%     TempLegend(i) = legend('Directed song', 'Undirected song');
%     set(TempLegend(i), 'Location', 'SouthEast', 'TextColor', 'w');
% end
% 
% %==========================================================================
% 
% %================ Plot match distribution with various potential thresholds
% %and spectrograms of matches near the threshold ===========================
% 
% figure;
% set(gcf, 'Color', 'k');
% set(gcf, 'Position', [128 216 1180 475]);
% Temp = cell2mat(Parameters.PreZeroThresholdUnDirResults{1});
% [SortedVals, SortedIndices] = sort(Temp(:,1));
% 
% Temp = Temp(SortedIndices, :);
% 
% XVals = linspace(0, max(Temp(:,1)), 100);
% plot(XVals, histc(Temp(:,1), XVals)/length(Temp(:,1)), 'w', 'LineWidth', 2);
% hold on;
% axis tight;
% TempAxis = axis;
% set(gca, 'Color', 'k');
% set(gca, 'XColor', 'w');
% set(gca, 'YColor', 'w');
% 
% % Plot threshold for other birds songs
% TempMatches = [];
% for i = 1:length(Parameters.OtherBirdSongResults),
%      TempMatches = [TempMatches; cell2mat(cellfun(@max, Parameters.OtherBirdSongResults{2}, 'UniformOutput', 0))];
% end
% Parameters.MaxOtherBirdSongs = max(TempMatches(:,1));
% plot(max(TempMatches(:,1))*ones(1,2), TempAxis(3:4), 'm--', 'LineWidth', 2);
% text(max(TempMatches(:,1))+0.05, 0.0003, ['Max match to other bird songs'], 'Rotation', 90, 'Color', 'm', 'FontSize', 14);
% 
% % Plot threshold for 100% shuffled songs
% plot(max(Parameters.MaxShuffledSongMatches)*ones(1,2), TempAxis(3:4), 'r--', 'LineWidth', 2);
% text(max(Parameters.MaxShuffledSongMatches)+0.05, 0.0003, ['Max match to 100% shuffled songs'], 'Rotation', 90, 'Color', 'r', 'FontSize', 14);
% 
% % First threshold for other birds songs
% for i = 1:length(Parameters.PartMaxShuffledSongMatches),
%     plot(max(Parameters.PartMaxShuffledSongMatches{i})*ones(1,2), TempAxis(3:4), 'b--', 'LineWidth', 2);
%     text(max(Parameters.PartMaxShuffledSongMatches{i})+0.05, 0.0003, ['Max match to ', num2str(25*i), '% shuffled songs'], 'Rotation', 90, 'Color', 'b', 'FontSize', 14);
% end
% annotation('textbox', [0.45 0.95 0.1 0.05], 'String', Parameters.BirdName, 'Color', 'w', 'FontSize', 20)
% axis([TempAxis(1:2) 0 0.001]);
% xlabel('Template match value for pre undir song', 'FontSize', 16, 'Color', 'w');
% set(gca, 'FontSize', 16);
% ylabel('%', 'FontSize', 16, 'Color', 'w');
% 
% figure;
% set(gcf, 'Color', 'k');
% set(gcf, 'Position', [141 6 900 685]);
% 
% % First plot spectrograms closest to other birds song max
% Index = find(Temp(:,1) > Parameters.MaxOtherBirdSongs, 1, 'first');
% 
% PlotIndex = 1;
% for i = Index-1:1:Index+1,
%     PlotSpectrogramInAxis([Parameters.PreDataDir{1}, FileSep], Parameters.PreUnDirSongFileNames{1}{Temp(i,3)}, Parameters.FileType, subplot(5, 3, PlotIndex), [(Temp(i,2)) (Temp(i,2) + 1)]);
%     PlotIndex = PlotIndex + 1;
%     title(Temp(i,1), 'Color', 'w', 'FontSize', 14);
% end
% annotation('textbox', [0.01 0.85 0.14 0.1], 'String', 'Matches near max of other bird songs', 'Color', 'm', 'FontSize', 16)
% 
% % First plot spectrograms closest to shuffled song max
% Index = find(Temp(:,1) > max(Parameters.MaxShuffledSongMatches), 1, 'first');
% 
% PlotIndex = 4;
% for i = Index-1:1:Index+1,
%     PlotSpectrogramInAxis([Parameters.PreDataDir{1}, FileSep], Parameters.PreUnDirSongFileNames{1}{Temp(i,3)}, Parameters.FileType, subplot(5, 3, PlotIndex), [(Temp(i,2)) (Temp(i,2) + 1)]);
%     PlotIndex = PlotIndex + 1;
%     title(Temp(i,1), 'Color', 'w', 'FontSize', 14);
% end
% annotation('textbox', [0.01 0.68 0.14 0.1], 'String', 'Matches near max of 100% shuffled songs', 'Color', 'r', 'FontSize', 16)
% 
% % First plot spectrograms closest to shuffled song max
% PlotIndex = 7;
% for j = length(Parameters.PartMaxShuffledSongMatches):-1:1,
%     Index = find(Temp(:,1) > max(Parameters.PartMaxShuffledSongMatches{j}), 1, 'first');
%     for i = Index-1:1:Index+1,
%         PlotSpectrogramInAxis([Parameters.PreDataDir{1}, FileSep], Parameters.PreUnDirSongFileNames{1}{Temp(i,3)}, Parameters.FileType, subplot(5, 3, PlotIndex), [(Temp(i,2)) (Temp(i,2) + 1)]);
%         PlotIndex = PlotIndex + 1;
%         title(Temp(i,1), 'Color', 'w', 'FontSize', 14);
%     end
%     annotation('textbox', [0.01 (0.68 - (4 - j)*0.17) 0.14 0.1], 'String', ['Matches near max of ', num2str(j*25), '% shuffled songs'], 'Color', 'b', 'FontSize', 16)
% end
% 
% 
% annotation('textbox', [0.45 0.95 0.1 0.05], 'String', Parameters.BirdName, 'Color', 'w', 'FontSize', 20)
% subplot(5, 3, 14);
% xlabel('Time - 1 sec window starting at the time of match', 'FontSize', 16, 'Color', 'w');
disp('Finished analysis');