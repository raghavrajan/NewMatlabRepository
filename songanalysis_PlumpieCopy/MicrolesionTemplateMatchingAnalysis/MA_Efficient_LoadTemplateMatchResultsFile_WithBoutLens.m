function [Results, ShuffleResults] = MA_Efficient_LoadTemplateMatchResultsFile_WithBoutLens(SongFile, OutputDir, TemplateFileName, ShufflePercentages, SongFileNo)

NumRepetition = 25;
Threshold = 0;

SmoothingKernelLen = 3;
SmoothingKernel = ones(SmoothingKernelLen,1)/SmoothingKernelLen;

[SyllableTemplateDir, SyllableTemplateFileName, SyllableTemplateExt] = fileparts(TemplateFileName);

SyllableTemplates = load(TemplateFileName);

for SyllTemp = 1:length(SyllableTemplates.SyllableTemplates),
    Labels{SyllTemp} = ['Syll_', SyllableTemplates.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Label];
end

TemplateMatchOutputDir = fullfile(OutputDir, [SyllableTemplateFileName, SyllableTemplateExt, '.TemplateMatchResults']);

for i = 1:length(Labels),
    % First load up the data for that file for the actual template matches
    SyllDir = fullfile(TemplateMatchOutputDir, Labels{i});
    if (~exist(fullfile(SyllDir, [SongFile, '.', Labels{i}, '.TempMatch.mat']), 'file'))
        Results = [];
        return;
    end

    Temp = load(fullfile(SyllDir, [SongFile, '.', Labels{i}, '.TempMatch.mat']));

    Results{i} = [];
    for j = 1:length(Temp.Bout),
        [Pks, Locs] = findpeaks(conv(Temp.Bout{j}.MaxBoutSeqMatch, SmoothingKernel, 'same'), 'MINPEAKHEIGHT', Threshold);
        Pks = Pks(:);
        Locs = Locs(:);
    
        Results{i} = [Results{i}; [Pks (Temp.Bout{j}.T(Locs)' + Temp.Bout{j}.BoutOnset/1000) ones(size(Pks))*SongFileNo]];
    end
    
    % Now load up the data from shuffled comparisons
    ShuffledComparisonDir = fullfile(TemplateMatchOutputDir, 'ShuffledSongComparisons');
    for j = 1:length(ShufflePercentages),
        ShuffleResults{i}{j} = [];
        for k = 1:NumRepetition,
            if (exist(fullfile(ShuffledComparisonDir, [SongFile, '.', Labels{i}, '.', num2str(ShufflePercentages(j)), '%.', num2str(k), '.TempMatch.mat']), 'file'))
                Temp = load(fullfile(ShuffledComparisonDir, [SongFile, '.', Labels{i}, '.', num2str(ShufflePercentages(j)), '%.', num2str(k), '.TempMatch.mat']));
                for BoutNo = 1:length(Temp.Bout),
                    ShuffleResults{i}{j}(end+1) = max(findpeaks(conv(Temp.Bout{BoutNo}.MaxBoutSeqMatch, SmoothingKernel, 'same')));
                end
            end
        end
    end
end       
disp('Finished');