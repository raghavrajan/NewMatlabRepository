function [Results] = MA_CompleteAnalysis(ParameterFile, SimilarityFigure, ConsistencyFigure)

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


%============= Load results of template matching ==========================
TemplateMatchOutputDir = '/home/raghav/MicrolesionAnalysisResults_MT_Templates/ShuffledSongComparisons/';
[MotifTemplateDir, MotifTemplateFileName, MotifTemplateExt] = fileparts(Parameters.MotifTemplateFileName);

TemplateMatchOutputDir = [TemplateMatchOutputDir, MotifTemplateFileName, MotifTemplateExt, '.TemplateMatchResults', FileSep];
disp('Loading results of shuffled song template matching ...');

NumberofFiles = min(50, length(Parameters.PreUnDirSongFileNames{1}));
RandomFiles = Parameters.PreUnDirSongFileNames{1}(1:NumberofFiles);

LabelCellArray = cellstr(char(ones(length(RandomFiles), 1)*double('Motif')));
OutputDirCellArray = cellstr(char(ones(length(RandomFiles), 1)*double(TemplateMatchOutputDir)));

Parameters.ShuffledSongMatches = cellfun(@MA_LoadShuffledSong_TemplateMatchResultsFile, RandomFiles, OutputDirCellArray, LabelCellArray, 'UniformOutput', 0);
Parameters.ShuffledSongMatches = Parameters.ShuffledSongMatches(:)';
Parameters.ShuffledSongMatches = cell2mat(Parameters.ShuffledSongMatches);
Parameters.ShuffledSongMatches = Parameters.ShuffledSongMatches(:);

Parameters.TemplateMatchThreshold = mean(Parameters.ShuffledSongMatches) + 3*std(Parameters.ShuffledSongMatches);

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
    Parameters.PartShuffledSongMatches{ShuffleIndex} = Parameters.PartShuffledSongMatches{ShuffleIndex}(:)';
    Parameters.PartShuffledSongMatches{ShuffleIndex} = cell2mat(Parameters.PartShuffledSongMatches{ShuffleIndex});
    Parameters.PartShuffledSongMatches{ShuffleIndex} = Parameters.PartShuffledSongMatches{ShuffleIndex}(:);
    Parameters.PartShuffleTemplateMatchThreshold{ShuffleIndex} = mean(Parameters.PartShuffledSongMatches{ShuffleIndex}) + 3*std(Parameters.PartShuffledSongMatches{ShuffleIndex});
end

% Choose one of the thresholds - 75% shuffled song 
ActualTemplateMatchThreshold = Parameters.PartShuffleTemplateMatchThreshold{3};

TemplateMatchOutputDir = '/home/raghav/MicrolesionAnalysisResults_MT_Templates/';
[MotifTemplateDir, MotifTemplateFileName, MotifTemplateExt] = fileparts(Parameters.MotifTemplateFileName);

TemplateMatchOutputDir = [TemplateMatchOutputDir, MotifTemplateFileName, MotifTemplateExt, '.TemplateMatchResults', FileSep];

disp('Loading results of song template matching ...');

% First for pre song
for i = 1:Parameters.NoPreDays,
    % First for dir song
    disp(['   Pre Day #', num2str(i), ' - directed song ...']); 
    LabelCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double('Motif')));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double(TemplateMatchOutputDir)));
    ThresholdCellArray = num2cell(ones(length(Parameters.PreDirSongFileNames{i}), 1) * ActualTemplateMatchThreshold);
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
    ThresholdCellArray = num2cell(ones(length(Parameters.PreUnDirSongFileNames{i}), 1) * ActualTemplateMatchThreshold);
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
    % First for dir song
    disp(['   Post Day #', num2str(i), ' - directed song ...']); 
    LabelCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double('Motif')));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double(TemplateMatchOutputDir)));
    ThresholdCellArray = num2cell(ones(length(Parameters.PostDirSongFileNames{i}), 1) * ActualTemplateMatchThreshold);
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
    ThresholdCellArray = num2cell(ones(length(Parameters.PostUnDirSongFileNames{i}), 1) * ActualTemplateMatchThreshold);
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
disp('Plotting match results ...');
figure;

for i = 1:length(Parameters.PreDirResults),
    DirAllMatches{i} = cell2mat(Parameters.PreDirResults{i});
    DirTotalSongTime{i} = cell2mat(Parameters.PreDirLens{i});
    
    if (~isempty(DirAllMatches{i}))
        subplot(2,1,1);
        errorbar(i, mean(DirAllMatches{i}(:,1)), std(DirAllMatches{i}(:,1)), 'ro');
        hold on;

        subplot(2,1,2);
        plot(i, size(DirAllMatches{i},1)*1000/sum(DirTotalSongTime{i}), 'ro');
        hold on;
    end
    
    UnDirAllMatches{i} = cell2mat(Parameters.PreUnDirResults{i});
    UnDirTotalSongTime{i} = cell2mat(Parameters.PreUnDirLens{i});
    
    if (~isempty(UnDirAllMatches{i}))
        subplot(2,1,1);
        errorbar(i, mean(UnDirAllMatches{i}(:,1)), std(UnDirAllMatches{i}(:,1)), 'bo');
        hold on;

        subplot(2,1,2);
        plot(i, size(UnDirAllMatches{i},1)*1000/sum(UnDirTotalSongTime{i}), 'bo');
        hold on;
    end
end
subplot(2,1,1);
hold on;
temp = axis;
plot([i i]+0.5, [0 temp(4)], 'k--', 'LineWidth', 2);

subplot(2,1,2);
hold on;
temp = axis;
plot([i i]+0.5, [0 temp(4)], 'k--', 'LineWidth', 2);

for i = 1:length(Parameters.PostDirResults),
    DirAllMatches{i + length(Parameters.PreDirResults)} = cell2mat(Parameters.PostDirResults{i});
    DirTotalSongTime{i + length(Parameters.PreDirResults)} = cell2mat(Parameters.PostDirLens{i});
    
    if (~isempty(DirAllMatches{i + length(Parameters.PreDirResults)}))
        subplot(2,1,1);
        errorbar(i + length(Parameters.PreDirResults), mean(DirAllMatches{i + length(Parameters.PreDirResults)}(:,1)), std(DirAllMatches{i + length(Parameters.PreDirResults)}(:,1)), 'ro');
        hold on;

        subplot(2,1,2);
        plot(i + length(Parameters.PreDirResults), size(DirAllMatches{i + length(Parameters.PreDirResults)},1)*1000/sum(DirTotalSongTime{i + length(Parameters.PreDirResults)}), 'ro');
        hold on;
    end
    
    UnDirAllMatches{i + length(Parameters.PreDirResults)} = cell2mat(Parameters.PostUnDirResults{i});
    UnDirTotalSongTime{i + length(Parameters.PreDirResults)} = cell2mat(Parameters.PostUnDirLens{i});
    
    if (~isempty(UnDirAllMatches{i + length(Parameters.PreDirResults)}))
        subplot(2,1,1);
        errorbar(i + length(Parameters.PreDirResults), mean(UnDirAllMatches{i + length(Parameters.PreDirResults)}(:,1)), std(UnDirAllMatches{i + length(Parameters.PreDirResults)}(:,1)), 'bo');
        hold on;

        subplot(2,1,2);
        plot(i + length(Parameters.PreDirResults), size(UnDirAllMatches{i + length(Parameters.PreDirResults)},1)*1000/sum(UnDirTotalSongTime{i + length(Parameters.PreDirResults)}), 'bo');
        hold on;
    end
end

subplot(2,1,1);
temp = axis;
plot([temp(1) temp(2)], [ActualTemplateMatchThreshold ActualTemplateMatchThreshold], 'k--', 'LineWidth', 2);

PostDay = inputdlg('Choose the post day that you want plotted', 'Post Day');
PostDay = str2double(PostDay{1}) + length(Parameters.PreDirResults);

figure(SimilarityFigure);
subplot(1,2,1);
scatter(mean(DirAllMatches{1}(:,1)), mean(UnDirAllMatches{1}(:,1)), Parameters.PercentTotalHVCremaining, 'k');
hold on;
plot([mean(DirAllMatches{1}(:,1)) mean(DirAllMatches{1}(:,1))], [(mean(UnDirAllMatches{1}(:,1)) + std(UnDirAllMatches{1}(:,1))) (mean(UnDirAllMatches{1}(:,1)) - std(UnDirAllMatches{1}(:,1)))], 'k');
plot([(mean(DirAllMatches{1}(:,1)) + std(DirAllMatches{1}(:,1))) (mean(DirAllMatches{1}(:,1)) - std(DirAllMatches{1}(:,1)))], [(mean(UnDirAllMatches{1}(:,1))) (mean(UnDirAllMatches{1}(:,1)))], 'k');

subplot(1,2,2);
scatter(mean(DirAllMatches{PostDay}(:,1)), mean(UnDirAllMatches{PostDay}(:,1)), Parameters.PercentTotalHVCremaining, 'k');
hold on;
plot([mean(DirAllMatches{PostDay}(:,1)) mean(DirAllMatches{PostDay}(:,1))], [(mean(UnDirAllMatches{PostDay}(:,1)) + std(UnDirAllMatches{PostDay}(:,1))) (mean(UnDirAllMatches{PostDay}(:,1)) - std(UnDirAllMatches{PostDay}(:,1)))], 'k');
plot([(mean(DirAllMatches{PostDay}(:,1)) + std(DirAllMatches{PostDay}(:,1))) (mean(DirAllMatches{PostDay}(:,1)) - std(DirAllMatches{PostDay}(:,1)))], [(mean(UnDirAllMatches{PostDay}(:,1))) (mean(UnDirAllMatches{PostDay}(:,1)))], 'k');

figure(ConsistencyFigure);
subplot(1,2,1);
scatter(size(DirAllMatches{1},1)*1000/sum(DirTotalSongTime{1}), size(UnDirAllMatches{1},1)*1000/sum(UnDirTotalSongTime{1}), Parameters.PercentTotalHVCremaining, 'k');
hold on;

subplot(1,2,2);
scatter(size(DirAllMatches{PostDay},1)*1000/sum(DirTotalSongTime{PostDay}), size(UnDirAllMatches{PostDay},1)*1000/sum(UnDirTotalSongTime{PostDay}), Parameters.PercentTotalHVCremaining, 'k');
hold on;

%==========================================================================

disp('Finished analysis');