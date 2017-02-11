function [Results] = MA_Analyze_TemplateMatchingTest(ParameterFile, NoteFileDir, RootOutputDir, SyllTemplateFile, RepetitionNum, ExampleNum)

%============== Some common variables =====================================
FileSep = filesep;
SmoothingKernelLen = 50;
SmoothingKernel = ones(SmoothingKernelLen,1)/SmoothingKernelLen;
MatchPrePostWindow = 0.005;
%==========================================================================

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

%============ Now to load up the results of syllable template matching =================
%============= Template Matching ==========================================

TemplateMatchOutputDir = [RootOutputDir, FileSep, 'PitchShiftedTemplates/'];

if (exist('SyllTemplateFile', 'var'))
    Parameters.SyllableTemplateFileName = SyllTemplateFile;
end

[SyllableTemplateDir, SyllableTemplateFileName, SyllableTemplateExt] = fileparts(Parameters.SyllableTemplateFileName);

TemplateMatchOutputDir = [TemplateMatchOutputDir, SyllableTemplateFileName, SyllableTemplateExt, '.TemplateMatchResults'];

ShuffledTemplateMatchOutputDir = fullfile(TemplateMatchOutputDir, 'ShuffledSongComparisons');

disp('Loading motif template ...');
Parameters.SyllableTemplate = load(Parameters.SyllableTemplateFileName);

for SyllTemp = 1:length(Parameters.SyllableTemplate.SyllableTemplates),
    Parameters.SyllableLabel{SyllTemp} = Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label;
    Parameters.SyllableTemplateLabel{SyllTemp} = ['Syll_', Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label];
    if (~exist([TemplateMatchOutputDir, FileSep, Parameters.SyllableTemplateLabel{SyllTemp}], 'dir'))
        mkdir([TemplateMatchOutputDir, FileSep, Parameters.SyllableTemplateLabel{SyllTemp}]);
    end
end

Parameters.ShufflePercentages = [25 100];

% Get all template Parameters
Results.Normalization = Parameters.SyllableTemplate.SyllableTemplates{1}{1}.MotifTemplate(1).Normalization;
Results.FFTWinSize = Parameters.SyllableTemplate.SyllableTemplates{1}{1}.MotifTemplate(1).FFTWinSize;
Results.FFTWinOverlap = Parameters.SyllableTemplate.SyllableTemplates{1}{1}.MotifTemplate(1).FFTWinOverlap;
Results.FreqIncrement = mean(diff(unique([Parameters.SyllableTemplate.SyllableTemplates{1}{1}.MotifTemplate.FreqStretch])));
Results.TimeIncrement = mean(diff(unique([Parameters.SyllableTemplate.SyllableTemplates{1}{1}.MotifTemplate.TimeStretch])));

% Now for each syllable, load up the file corresponding to the template
% match values

% All of this testing has been done with undirected song, so I can do this
% for the 1st day of undirected singing.
% Also, did all the shuffled matching only for the first song bout, so will
% consider only the first song bout now.

for i = 1:length(Parameters.PreUnDirSongFileNames{1}),
    % First load up the note file with info about the labelled syllables
    NoteFileInfo = load(fullfile(NoteFileDir, [Parameters.PreUnDirSongFileNames{1}{i}, '.not.mat']));
    % First load up all the syllable matches
    for j = 1:length(Parameters.SyllableTemplateLabel),
        TemplateMatchResultsFile = [Parameters.PreUnDirSongFileNames{1}{i}, '.', Parameters.SyllableTemplateLabel{j}, '.TempMatch.mat'];
        TemplateMatchResults = load(fullfile(TemplateMatchOutputDir, fullfile(Parameters.SyllableTemplateLabel{j}, TemplateMatchResultsFile)));
        
        xx{j} = TemplateMatchResults.Bout{1}.T(1):0.0001:TemplateMatchResults.Bout{1}.T(length(TemplateMatchResults.Bout{1}.MaxBoutSeqMatch));
        SmoothedResults = spline(TemplateMatchResults.Bout{1}.T(1:length(TemplateMatchResults.Bout{1}.MaxBoutSeqMatch)), TemplateMatchResults.Bout{1}.MaxBoutSeqMatch, xx{j});
        SmoothedTemplateMatchResults{j} = conv(SmoothedResults, SmoothingKernel, 'same');
        xx{j} = xx{j} + TemplateMatchResults.Bout{1}.BoutOnset/1000;
    end
    
    % Now to find best syllable and best other syllable matches at location
    % of labelled syllable of each type.
    for j = 1:length(Parameters.SyllableTemplateLabel),
        Matches = find(NoteFileInfo.labels == Parameters.SyllableLabel{j});
        SyllableMatches = [];
        OtherSyllableMatches = [];
        OtherSyllableMatchIdentity = [];
        for k = 1:length(Parameters.SyllableTemplateLabel),
            for MatchNo = 1:length(Matches),
                CorrespondingMatchLocations = find((xx{k} >= (NoteFileInfo.onsets(Matches(MatchNo))/1000 - MatchPrePostWindow)) & (xx{k} <= (NoteFileInfo.onsets(Matches(MatchNo))/1000 + MatchPrePostWindow)));
                if (j == k)
                    SyllableMatches = [SyllableMatches; max(SmoothedTemplateMatchResults{k}(CorrespondingMatchLocations))];
                else
                    OtherSyllableMatches = [OtherSyllableMatches; max(SmoothedTemplateMatchResults{k}(CorrespondingMatchLocations))];
                    OtherSyllableMatchIdentity = [OtherSyllableMatchIdentity; k];
                end
            end
        end
        Results.OtherSyllableMatches{j} = OtherSyllableMatches;
        Results.OtherSyllableMatchIdentities{j} = OtherSyllableMatchIdentity;
        Results.SyllableMatches{j} = SyllableMatches;
    end
    
    % Now to get results of shuffled matching - just get the maximum value
    % for each of the files
    for j = 1:length(Parameters.SyllableTemplateLabel),
        for k = 1:length(Parameters.ShufflePercentages),
            for Repetition = 1:RepetitionNum,
                TemplateMatchResultsFile = [Parameters.PreUnDirSongFileNames{1}{i}, '.', Parameters.SyllableTemplateLabel{j}, '.', num2str(Parameters.ShufflePercentages(k)), '%.', num2str(Repetition), '.TempMatch.mat'];
                TemplateMatchResults = load(fullfile(ShuffledTemplateMatchOutputDir, TemplateMatchResultsFile));
                Temp_xx = TemplateMatchResults.Bout{1}.T(1):0.0001:TemplateMatchResults.Bout{1}.T(length(TemplateMatchResults.Bout{1}.MaxBoutSeqMatch));
                SmoothedResults = spline(TemplateMatchResults.Bout{1}.T(1:length(TemplateMatchResults.Bout{1}.MaxBoutSeqMatch)), TemplateMatchResults.Bout{1}.MaxBoutSeqMatch, Temp_xx);
                Temp_SmoothedTemplateMatchResults = conv(SmoothedResults, SmoothingKernel, 'same');
                Results.ShuffledTemplateMatchResults{j}{k}(Repetition) = max(Temp_SmoothedTemplateMatchResults);
            end
        end
    end
end

%==========================================================================

disp('Finished analysis');
