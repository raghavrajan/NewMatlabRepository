function [] = MA_CompleteSyllableAnalysis_EfficientVersion(ParameterFile, NoteFileDir, RootOutputDir, varargin)

%==========================================================================
% Function for analysing the effects of a treatment on song - written in
% the context of the analysis of effects of HVC microlesions. Specifically
% for analysing individual syllables
% Raghav Rajan - 01st July 2014
%==========================================================================

if (nargin > 3)
    SyllTemplateFile = varargin{1};
end

%============== Some common variables =====================================
FileSep = filesep;
OutputDir = fullfile(RootOutputDir, 'MicrolesionAnalysisResults');

if (~exist(OutputDir, 'dir'))
    mkdir(RootOutputDir, 'MicrolesionAnalysisResults');
end

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

%======Now load up note files with bout information =======================
% Load up note files which have bouts labelled with a 'B'. This has been
% done to get a better idea of consistency of song production.

disp('Loading up note files ...');

if (~isempty(NoteFileDir))
    % First for pre days
    for i = 1:Parameters.NoPreDays,
        disp(['   Pre Day #', num2str(i), ' - directed song ...']); 
        % First directed songs
        NoteFileDirCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double(NoteFileDir)));

        % using cellfun so that i iterate over each element of the cell array.
        % To use cellfun, all of the other inputs also have to be in the form
        % of cell arrays of the same length - so the previous three lines
        % convert file type, data dir and output dir - common parameters for
        % all of the files into cell arrays

        [Parameters.PreDirBoutOnsets{i}, Parameters.PreDirBoutOffsets{i}, Parameters.PreDirBoutLens{i}] = cellfun(@MA_LoadBoutNoteFiles, Parameters.PreDirSongFileNames{i}, NoteFileDirCellArray, 'UniformOutput', 0);

        % next undirected songs
        disp(['   Pre Day #', num2str(i), ' - undirected song ...']); 
        NoteFileDirCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double(NoteFileDir)));
        [Parameters.PreUnDirBoutOnsets{i}, Parameters.PreUnDirBoutOffsets{i}, Parameters.PreUnDirBoutLens{i}] = cellfun(@MA_LoadBoutNoteFiles, Parameters.PreUnDirSongFileNames{i}, NoteFileDirCellArray, 'UniformOutput', 0);
    end

    % Next for post days
    for i = 1:Parameters.NoPostDays,
        % First directed songs
        disp(['   Post Day #', num2str(i), ' - directed song ...']); 
        
        SlashIndex = find(Parameters.PostDataDir{i} == FileSep);
        NoteFileDirCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double(NoteFileDir)));

        % using cellfun so that i iterate over each element of the cell array.
        % To use cellfun, all of the other inputs also have to be in the form
        % of cell arrays of the same length - so the previous three lines
        % convert file type, data dir and output dir - common parameters for
        % all of the files into cell arrays

        [Parameters.PostDirBoutOnsets{i}, Parameters.PostDirBoutOffsets{i}, Parameters.PostDirBoutLens{i}] = cellfun(@MA_LoadBoutNoteFiles, Parameters.PostDirSongFileNames{i}, NoteFileDirCellArray, 'UniformOutput', 0);

        % next undirected songs
        disp(['   Post Day #', num2str(i), ' - undirected song ...']); 
        NoteFileDirCellArray = cellstr(char(ones(length(Parameters.PostUnDirSongFileNames{i}), 1)*double(NoteFileDir)));
        [Parameters.PostUnDirBoutOnsets{i}, Parameters.PostUnDirBoutOffsets{i}, Parameters.PostUnDirBoutLens{i}] = cellfun(@MA_LoadBoutNoteFiles, Parameters.PostUnDirSongFileNames{i}, NoteFileDirCellArray, 'UniformOutput', 0);
    end
end
%==========================================================================

%======Now segment files (Aronov Fee style) ===============================
% Segment each file separately, minimum interval = 7ms, minimum duration =
% 7ms and the thresholds are determined separately for each file similar to
% Aronov and Fee (J.Neurosci) paper.
% Use the bout information from the previous thing to segment only the bout
% part of the song

disp('Segmenting files into note files ...');

% First for pre days
for i = 1:Parameters.NoPreDays,
    disp(['   Pre Day #', num2str(i), ' - directed song ...']); 
    % First directed songs
    FileTypeCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double(Parameters.FileType)));
    RawDataDirCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double(Parameters.PreDataDir{i})));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double([OutputDir, Parameters.BirdName, '.NoteFiles', FileSep])));
    BoutCellArray = num2cell(ones(length(Parameters.PreDirSongFileNames{i}), 1));
    % using cellfun so that i iterate over each element of the cell array.
    % To use cellfun, all of the other inputs also have to be in the form
    % of cell arrays of the same length - so the previous three lines
    % convert file type, data dir and output dir - common parameters for
    % all of the files into cell arrays
    
    [Parameters.PreDirOnsets{i}, Parameters.PreDirOffsets{i}, Parameters.PreDirSyllDurs{i}, Parameters.PreDirGapDurs{i}, Parameters.PreDirThresholds{i}, Parameters.PreDirLens{i}, Parameters.PreDirFs{i}] = cellfun(@MA_SegmentFiles, Parameters.PreDirSongFileNames{i}, FileTypeCellArray, RawDataDirCellArray, OutputDirCellArray, Parameters.PreDirBoutOnsets{i}, Parameters.PreDirBoutOffsets{i}, BoutCellArray, 'UniformOutput', 0);

    % next undirected songs
    disp(['   Pre Day #', num2str(i), ' - undirected song ...']); 
    FileTypeCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double(Parameters.FileType)));
    RawDataDirCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double(Parameters.PreDataDir{i})));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double([OutputDir, Parameters.BirdName, '.NoteFiles', FileSep])));
    BoutCellArray = num2cell(ones(length(Parameters.PreUnDirSongFileNames{i}), 1));
    [Parameters.PreUnDirOnsets{i}, Parameters.PreUnDirOffsets{i}, Parameters.PreUnDirSyllDurs{i}, Parameters.PreUnDirGapDurs{i}, Parameters.PreUnDirThresholds{i}, Parameters.PreUnDirLens{i}, Parameters.PreUnDirFs{i}] = cellfun(@MA_SegmentFiles, Parameters.PreUnDirSongFileNames{i}, FileTypeCellArray, RawDataDirCellArray, OutputDirCellArray, Parameters.PreUnDirBoutOnsets{i}, Parameters.PreUnDirBoutOffsets{i}, BoutCellArray, 'UniformOutput', 0);
end

% Next for post days
for i = 1:Parameters.NoPostDays,
    % First directed songs
    disp(['   Post Day #', num2str(i), ' - directed song ...']); 
    FileTypeCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double(Parameters.FileType)));
    RawDataDirCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double(Parameters.PostDataDir{i})));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double([OutputDir, Parameters.BirdName, '.NoteFiles', FileSep])));
    BoutCellArray = num2cell(ones(length(Parameters.PostDirSongFileNames{i}), 1));
    
    % using cellfun so that i iterate over each element of the cell array.
    % To use cellfun, all of the other inputs also have to be in the form
    % of cell arrays of the same length - so the previous three lines
    % convert file type, data dir and output dir - common parameters for
    % all of the files into cell arrays
    
    [Parameters.PostDirOnsets{i}, Parameters.PostDirOffsets{i}, Parameters.PostDirSyllDurs{i}, Parameters.PostDirGapDurs{i}, Parameters.PostDirThresholds{i}, Parameters.PostDirLens{i}, Parameters.PostDirFs{i}] = cellfun(@MA_SegmentFiles, Parameters.PostDirSongFileNames{i}, FileTypeCellArray, RawDataDirCellArray, OutputDirCellArray, Parameters.PostDirBoutOnsets{i}, Parameters.PostDirBoutOffsets{i}, BoutCellArray, 'UniformOutput', 0);

    % next undirected songs
    disp(['   Post Day #', num2str(i), ' - undirected song ...']); 
    FileTypeCellArray = cellstr(char(ones(length(Parameters.PostUnDirSongFileNames{i}), 1)*double(Parameters.FileType)));
    RawDataDirCellArray = cellstr(char(ones(length(Parameters.PostUnDirSongFileNames{i}), 1)*double(Parameters.PostDataDir{i})));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PostUnDirSongFileNames{i}), 1)*double([OutputDir, Parameters.BirdName, '.NoteFiles', FileSep])));
    BoutCellArray = num2cell(ones(length(Parameters.PostUnDirSongFileNames{i}), 1));
    [Parameters.PostUnDirOnsets{i}, Parameters.PostUnDirOffsets{i}, Parameters.PostUnDirSyllDurs{i}, Parameters.PostUnDirGapDurs{i}, Parameters.PostUnDirThresholds{i}, Parameters.PostUnDirLens{i}, Parameters.PostUnDirFs{i}] = cellfun(@MA_SegmentFiles, Parameters.PostUnDirSongFileNames{i}, FileTypeCellArray, RawDataDirCellArray, OutputDirCellArray, Parameters.PostUnDirBoutOnsets{i}, Parameters.PostUnDirBoutOffsets{i}, BoutCellArray, 'UniformOutput', 0);
end
%==========================================================================

%======Now plot syllable and gap duration histograms ======================
% disp('Plotting syllable and gap duration histograms ...');
% Colours = 'rgbcmy';
% Symbols = 'o*';
% 
% BinSize = 5;
% MaxSyllDur = 350;
% MaxGapDur = 150;
% 
% SyllEdges = 0:BinSize:MaxSyllDur;
% GapEdges = 0:BinSize:MaxGapDur;
% 
% SyllGapDurHistFigure = figure;
% 
% Legend = [];
% % First for pre days
% for i = 1:Parameters.NoPreDays,
%     Legend{end+1} = ['Pre #', num2str(i)];
%     % First for directed song syllable durations
%     subplot(2, 2, 1);
%     hold on;
%     NonZeroDurs = find(cellfun(@length, Parameters.PreDirSyllDurs{i}));
%     plot(SyllEdges, histc(cell2mat(Parameters.PreDirSyllDurs{i}(NonZeroDurs)), SyllEdges) * 100 / sum(histc(cell2mat(Parameters.PreDirSyllDurs{i}(NonZeroDurs)), SyllEdges)), [Colours(i), Symbols(1), '-']);
%     
%     % Next for directed song gap durations
%     subplot(2, 2, 2);
%     hold on;
%     NonZeroDurs = find(cellfun(@length, Parameters.PreDirGapDurs{i}));
%     plot(GapEdges, histc(cell2mat(Parameters.PreDirGapDurs{i}(NonZeroDurs)), GapEdges) * 100 / sum(histc(cell2mat(Parameters.PreDirGapDurs{i}(NonZeroDurs)), GapEdges)), [Colours(i), Symbols(1), '-']);
% 
%     % Next for undirected song syllable durations
%     subplot(2, 2, 3);
%     hold on;
%     NonZeroDurs = find(cellfun(@length, Parameters.PreUnDirSyllDurs{i}));
%     plot(SyllEdges, histc(cell2mat(Parameters.PreUnDirSyllDurs{i}(NonZeroDurs)), SyllEdges) * 100 / sum(histc(cell2mat(Parameters.PreUnDirSyllDurs{i}(NonZeroDurs)), SyllEdges)), [Colours(i), Symbols(1), '-']);
%     
%     % Next for undirected song gap durations
%     subplot(2, 2, 4);
%     hold on;
%     NonZeroDurs = find(cellfun(@length, Parameters.PreUnDirGapDurs{i}));
%     plot(GapEdges, histc(cell2mat(Parameters.PreUnDirGapDurs{i}(NonZeroDurs)), GapEdges) * 100 / sum(histc(cell2mat(Parameters.PreUnDirGapDurs{i}(NonZeroDurs)), GapEdges)), [Colours(i), Symbols(1), '-']);
% end
% 
% % Next for post days
% for i = 1:Parameters.NoPostDays,
%     Legend{end+1} = ['Post #', num2str(i)];
%     % First for directed song syllable durations
%     subplot(2, 2, 1);
%     hold on;
%     NonZeroDurs = find(cellfun(@length, Parameters.PostDirSyllDurs{i}));
%     plot(SyllEdges, histc(cell2mat(Parameters.PostDirSyllDurs{i}(NonZeroDurs)), SyllEdges) * 100 / sum(histc(cell2mat(Parameters.PostDirSyllDurs{i}(NonZeroDurs)), SyllEdges)), [Colours(i), Symbols(2), '-']);
%     
%     % Next for directed song gap durations
%     subplot(2, 2, 2);
%     hold on;
%     NonZeroDurs = find(cellfun(@length, Parameters.PostDirGapDurs{i}));
%     plot(GapEdges, histc(cell2mat(Parameters.PostDirGapDurs{i}(NonZeroDurs)), GapEdges) * 100 / sum(histc(cell2mat(Parameters.PostDirGapDurs{i}(NonZeroDurs)), GapEdges)), [Colours(i), Symbols(2), '-']);
% 
%     % Next for undirected song syllable durations
%     subplot(2, 2, 3);
%     hold on;
%     NonZeroDurs = find(cellfun(@length, Parameters.PostUnDirSyllDurs{i}));
%     plot(SyllEdges, histc(cell2mat(Parameters.PostUnDirSyllDurs{i}(NonZeroDurs)), SyllEdges) * 100 / sum(histc(cell2mat(Parameters.PostUnDirSyllDurs{i}(NonZeroDurs)), SyllEdges)), [Colours(i), Symbols(2), '-']);
%     
%     % Next for undirected song gap durations
%     subplot(2, 2, 4);
%     hold on;
%     NonZeroDurs = find(cellfun(@length, Parameters.PostUnDirGapDurs{i}));
%     plot(GapEdges, histc(cell2mat(Parameters.PostUnDirGapDurs{i}(NonZeroDurs)), GapEdges) * 100 / sum(histc(cell2mat(Parameters.PostUnDirGapDurs{i}(NonZeroDurs)), GapEdges)), [Colours(i), Symbols(2), '-']);
% end
% 
% figure(SyllGapDurHistFigure);
% set(gcf, 'Color', 'k', 'Position', [400 150 1000 600]);
% annotation('textbox', [0.45 0.95 0.1 0.05], 'String', Parameters.BirdName, 'Color', 'w', 'FontSize', 20)
% 
% for i = 1:4,
%     subplot(2, 2, i);
%     set(gca, 'Color', 'k');
%     set(gca, 'XColor', 'w');
%     set(gca, 'YColor', 'w');
%     axis tight;
%     Axis(i,:) = axis;
% end
% Axis(:,[1 3]) = 0;
% Axis([1 3], 2) = MaxSyllDur;
% Axis([2 4], 2) = MaxGapDur;
% Axis([1 3], 4) = max(Axis([1 3], 4));
% Axis([2 4], 4) = max(Axis([2 4], 4));
% 
% Title{1} = 'Dir - syllable durations';
% Title{2} = 'Dir - gap durations';
% Title{3} = 'Undir - syllable durations';
% Title{4} = 'Undir - gap durations';
% 
% for i = 1:4,
%     subplot(2, 2, i);
%     legend(Legend, 'TextColor', 'w', 'FontSize', 14);
%     set(gca, 'Box', 'on');
%     axis(Axis(i,:));
%     if (i > 2)
%         xlabel('Duration (ms)', 'FontSize', 14);
%     end
%     
%     if (mod(i, 2) == 1)
%         ylabel('%', 'FontSize', 14);
%     end
%     title(Title{i}, 'FontSize', 14, 'Color', 'w');
% end
%==========================================================================

%============ Now to start the syllable template matching =================
%============= Template Matching ==========================================

TemplateMatchOutputDir = [RootOutputDir, FileSep, 'PitchShiftedTemplates/'];

if (~exist(TemplateMatchOutputDir, 'dir'))
    mkdir(RootOutputDir, 'PitchShiftedTemplates');
end

if (exist('SyllTemplateFile', 'var'))
    Parameters.SyllableTemplateFileName = SyllTemplateFile;
end

[SyllableTemplateDir, SyllableTemplateFileName, SyllableTemplateExt] = fileparts(Parameters.SyllableTemplateFileName);

TemplateMatchOutputDir = [TemplateMatchOutputDir, SyllableTemplateFileName, SyllableTemplateExt, '.TemplateMatchResults'];
if (~exist(TemplateMatchOutputDir, 'dir'))
    mkdir(TemplateMatchOutputDir);
end

if (~exist(fullfile(TemplateMatchOutputDir, 'ShuffledSongComparisons'), 'dir'))
    mkdir(fullfile(TemplateMatchOutputDir, 'ShuffledSongComparisons'));
end

disp('Loading motif template ...');
Parameters.SyllableTemplate = load(Parameters.SyllableTemplateFileName);

for SyllTemp = 1:length(Parameters.SyllableTemplate.SyllableTemplates),
    Parameters.SyllableTemplateLabel{SyllTemp} = ['Syll_', Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label];
    if (~exist([TemplateMatchOutputDir, FileSep, Parameters.SyllableTemplateLabel{SyllTemp}], 'dir'))
        mkdir([TemplateMatchOutputDir, FileSep, Parameters.SyllableTemplateLabel{SyllTemp}]);
    end
end

disp('Doing template matching ...');

Parameters.ShufflePercentages = [25 100];

% First for pre song
for i = 1:Parameters.NoPreDays,
    % First for dir song
    disp(['      Pre Day #', num2str(i), ' - directed song ...']); 
    
    FileTypeCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double(Parameters.FileType)));
    RawDataDirCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double(Parameters.PreDataDir{i})));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double([TemplateMatchOutputDir, FileSep])));
    TemplateTypeCellArray = cellstr(char(ones(length(Parameters.PreDirSongFileNames{i}), 1)*double('Spectrogram')));
    
    LabelCellArray = cell(size(Parameters.PreDirSongFileNames{i}));
    SyllableTemplateCellArray = cell(size(Parameters.PreDirSongFileNames{i}));
    ShufflePercentageCellArray = cell(size(Parameters.PreDirSongFileNames{i}));

    FileIndices = 1:1:length(Parameters.PreDirSongFileNames{i});
    FileIndices = FileIndices(:);
    FileIndexCellArray = num2cell(FileIndices);
    
    for j = 1:length(SyllableTemplateCellArray),
        LabelCellArray{j} = Parameters.SyllableTemplateLabel;
        SyllableTemplateCellArray{j} = Parameters.SyllableTemplate.SyllableTemplates;
        ShufflePercentageCellArray{j} = Parameters.ShufflePercentages;
    end
    
    % using cellfun so that i iterate over each element of the cell array.
    % To use cellfun, all of the other inputs also have to be in the form
    % of cell arrays of the same length - so the previous three lines
    % convert file type, data dir and output dir - common parameters for
    % all of the files into cell arrays
    cellfun(@MA_Efficient_TemplateMatch_WithBoutLens, RawDataDirCellArray, Parameters.PreDirSongFileNames{i}, FileTypeCellArray, SyllableTemplateCellArray, LabelCellArray, ShufflePercentageCellArray, OutputDirCellArray, TemplateTypeCellArray, Parameters.PreDirBoutOnsets{i}, Parameters.PreDirBoutOffsets{i}, FileIndexCellArray, 'UniformOutput', 0);

    % next undirected songs
    disp(['      Pre Day #', num2str(i), ' - undirected song ...']); 

    FileTypeCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double(Parameters.FileType)));
    RawDataDirCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double(Parameters.PreDataDir{i})));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double([TemplateMatchOutputDir, FileSep])));
    TemplateTypeCellArray = cellstr(char(ones(length(Parameters.PreUnDirSongFileNames{i}), 1)*double('Spectrogram')));
    
    LabelCellArray = cell(size(Parameters.PreUnDirSongFileNames{i}));
    SyllableTemplateCellArray = cell(size(Parameters.PreUnDirSongFileNames{i}));
    ShufflePercentageCellArray = cell(size(Parameters.PreUnDirSongFileNames{i}));

    FileIndices = 1:1:length(Parameters.PreUnDirSongFileNames{i});
    FileIndices = FileIndices(:);
    FileIndexCellArray = num2cell(FileIndices);

    for j = 1:length(SyllableTemplateCellArray),
        LabelCellArray{j} = Parameters.SyllableTemplateLabel;
        SyllableTemplateCellArray{j} = Parameters.SyllableTemplate.SyllableTemplates;
        ShufflePercentageCellArray{j} = Parameters.ShufflePercentages;
    end
    % using cellfun so that i iterate over each element of the cell array.
    % To use cellfun, all of the other inputs also have to be in the form
    % of cell arrays of the same length - so the previous three lines
    % convert file type, data dir and output dir - common parameters for
    % all of the files into cell arrays
    cellfun(@MA_Efficient_TemplateMatch_WithBoutLens, RawDataDirCellArray, Parameters.PreUnDirSongFileNames{i}, FileTypeCellArray, SyllableTemplateCellArray, LabelCellArray, ShufflePercentageCellArray, OutputDirCellArray, TemplateTypeCellArray, Parameters.PreUnDirBoutOnsets{i}, Parameters.PreUnDirBoutOffsets{i}, FileIndexCellArray, 'UniformOutput', 0);
end

% Next for post song
for i = 1:Parameters.NoPostDays,
    % First for dir song
    disp(['      Post Day #', num2str(i), ' - directed song ...']); 
    
    FileTypeCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double(Parameters.FileType)));
    RawDataDirCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double(Parameters.PostDataDir{i})));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double([TemplateMatchOutputDir, FileSep])));
    TemplateTypeCellArray = cellstr(char(ones(length(Parameters.PostDirSongFileNames{i}), 1)*double('Spectrogram')));
    
    LabelCellArray = cell(size(Parameters.PostDirSongFileNames{i}));
    SyllableTemplateCellArray = cell(size(Parameters.PostDirSongFileNames{i}));
    ShufflePercentageCellArray = cell(size(Parameters.PostDirSongFileNames{i}));

    FileIndices = 1:1:length(Parameters.PostDirSongFileNames{i});
    FileIndices = FileIndices(:);
    FileIndexCellArray = num2cell(FileIndices);
    
    for j = 1:length(SyllableTemplateCellArray),
        LabelCellArray{j} = Parameters.SyllableTemplateLabel;
        SyllableTemplateCellArray{j} = Parameters.SyllableTemplate.SyllableTemplates;
        ShufflePercentageCellArray{j} = Parameters.ShufflePercentages;
    end
    

    % using cellfun so that i iterate over each element of the cell array.
    % To use cellfun, all of the other inputs also have to be in the form
    % of cell arrays of the same length - so the previous three lines
    % convert file type, data dir and output dir - common parameters for
    % all of the files into cell arrays
    cellfun(@MA_Efficient_TemplateMatch_WithBoutLens, RawDataDirCellArray, Parameters.PostDirSongFileNames{i}, FileTypeCellArray, SyllableTemplateCellArray, LabelCellArray, ShufflePercentageCellArray, OutputDirCellArray, TemplateTypeCellArray, Parameters.PostDirBoutOnsets{i}, Parameters.PostDirBoutOffsets{i}, FileIndexCellArray, 'UniformOutput', 0);

    % next undirected songs
    disp(['      Post Day #', num2str(i), ' - undirected song ...']); 

    FileTypeCellArray = cellstr(char(ones(length(Parameters.PostUnDirSongFileNames{i}), 1)*double(Parameters.FileType)));
    RawDataDirCellArray = cellstr(char(ones(length(Parameters.PostUnDirSongFileNames{i}), 1)*double(Parameters.PostDataDir{i})));
    OutputDirCellArray = cellstr(char(ones(length(Parameters.PostUnDirSongFileNames{i}), 1)*double([TemplateMatchOutputDir, FileSep])));
    TemplateTypeCellArray = cellstr(char(ones(length(Parameters.PostUnDirSongFileNames{i}), 1)*double('Spectrogram')));
    
    LabelCellArray = cell(size(Parameters.PostUnDirSongFileNames{i}));
    SyllableTemplateCellArray = cell(size(Parameters.PostUnDirSongFileNames{i}));
    ShufflePercentageCellArray = cell(size(Parameters.PostUnDirSongFileNames{i}));

    FileIndices = 1:1:length(Parameters.PostUnDirSongFileNames{i});
    FileIndices = FileIndices(:);
    FileIndexCellArray = num2cell(FileIndices);

    for j = 1:length(SyllableTemplateCellArray),
        LabelCellArray{j} = Parameters.SyllableTemplateLabel;
        SyllableTemplateCellArray{j} = Parameters.SyllableTemplate.SyllableTemplates;
        ShufflePercentageCellArray{j} = Parameters.ShufflePercentages;
    end
    % using cellfun so that i iterate over each element of the cell array.
    % To use cellfun, all of the other inputs also have to be in the form
    % of cell arrays of the same length - so the previous three lines
    % convert file type, data dir and output dir - common parameters for
    % all of the files into cell arrays
    cellfun(@MA_Efficient_TemplateMatch_WithBoutLens, RawDataDirCellArray, Parameters.PostUnDirSongFileNames{i}, FileTypeCellArray, SyllableTemplateCellArray, LabelCellArray, ShufflePercentageCellArray, OutputDirCellArray, TemplateTypeCellArray, Parameters.PostUnDirBoutOnsets{i}, Parameters.PostUnDirBoutOffsets{i}, FileIndexCellArray, 'UniformOutput', 0);
end

%==========================================================================

disp('Finished analysis');