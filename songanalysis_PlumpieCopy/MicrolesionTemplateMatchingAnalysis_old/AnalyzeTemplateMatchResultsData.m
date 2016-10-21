function [Results] = AnalyzeTemplateMatchResultsData(handles)

SmoothingKernelLen = 1;
SmoothingKernel = ones(SmoothingKernelLen,1)/SmoothingKernelLen;

MPH = 0;
for i = 1:handles.NoofDaysToAnalyse,
    % First analyse the output data from directed songs

    Fid = fopen(handles.DirSongFileList{i}, 'r');
    Temp = textscan(Fid, '%s', 'delimiter', '\n');
    fclose(Fid);

    SongFiles = Temp{1};

    MotifPeaks = [];
    DirFileLens = [];
    for j = 1:length(handles.SyllTemplates.SyllableTemplates),
        SyllPeaks{j} = [];
    end
    
    for SongFileNo = 1:length(SongFiles),
        SongFile = SongFiles{SongFileNo};
        Slash = find((SongFile == '/') | (SongFile == '\'));
        if (~isempty(Slash))
            SongFile = SongFile(Slash(end)+1:end);
        end
        cd(handles.DirNormal_MotifTemplateMatchOutputFilesDir);
        Temp = load([SongFile, '.Motif.TempMatch.mat']);
        DirFileLens = [DirFileLens; Temp.Bout.FileLength];
        
        [Pks, Locs] = findpeaks(conv(Temp.Bout.MaxBoutSeqMatch, SmoothingKernel, 'same'), 'MINPEAKHEIGHT', MPH);
        Pks = Pks(:);
        Locs = Locs(:);
        
        MotifPeaks = [MotifPeaks; [Pks Temp.Bout.T(Locs)' ones(size(Pks))*SongFileNo]];
        
        cd(handles.DirNormal_SyllTemplateMatchOutputFilesDir);
        for j = 1:length(handles.SyllTemplates.SyllableTemplates),
            Temp = load([SongFile, '.', handles.SyllTemplates.SyllableTemplates{j}{1}.MotifTemplate(1).Label, '.TempMatch.mat']);
            [Pks, Locs] = findpeaks(conv(Temp.Bout.MaxBoutSeqMatch, SmoothingKernel, 'same'), 'MINPEAKHEIGHT', MPH);
            Pks = Pks(:);
            Locs = Locs(:);
        
            SyllPeaks{j} = [SyllPeaks{j}; [Pks Temp.Bout.T(Locs)' ones(size(Pks))*SongFileNo]];
        end
    end
    
    Results(i).DirSongFiles = SongFiles;
    Results(i).DirFileLens = DirFileLens;
    
    Results(i).DirPeaks{1} = MotifPeaks;
    Results(i).Labels{1} = handles.MotifTemplate.MotifTemplate(1).Label;
    
    if (exist('SyllPeaks', 'var'))
        for j = 1:length(SyllPeaks),
            Results(i).DirPeaks{j+1} = SyllPeaks{j};
            Results(i).Labels{j+1} = handles.SyllTemplates.SyllableTemplates{j}{1}.MotifTemplate(1).Label;
        end
    end
    
    % Next analyse the output data from undirected songs
    Fid = fopen(handles.UnDirSongFileList{i}, 'r');
    Temp = textscan(Fid, '%s', 'delimiter', '\n');
    fclose(Fid);

    SongFiles = Temp{1};

    UnDirFileLens = [];
    MotifPeaks = [];
    for j = 1:length(handles.SyllTemplates.SyllableTemplates),
        SyllPeaks{j} = [];
    end
    
    for SongFileNo = 1:length(SongFiles),
        SongFile = SongFiles{SongFileNo};
        Slash = find((SongFile == '/') | (SongFile == '\'));
        if (~isempty(Slash))
            SongFile = SongFile(Slash(end)+1:end);
        end
        cd(handles.UnDirNormal_MotifTemplateMatchOutputFilesDir);
        Temp = load([SongFile, '.Motif.TempMatch.mat']);
        
        UnDirFileLens = [UnDirFileLens; Temp.Bout.FileLength];
        [Pks, Locs] = findpeaks(conv(Temp.Bout.MaxBoutSeqMatch, SmoothingKernel, 'same'), 'MINPEAKHEIGHT', MPH);
        Pks = Pks(:);
        Locs = Locs(:);
        
        MotifPeaks = [MotifPeaks; [Pks Temp.Bout.T(Locs)' ones(size(Pks))*SongFileNo]];
        
        cd(handles.UnDirNormal_SyllTemplateMatchOutputFilesDir);
        clear TempBout;
        for j = 1:length(handles.SyllTemplates.SyllableTemplates),
            Temp = load([SongFile, '.', handles.SyllTemplates.SyllableTemplates{j}{1}.MotifTemplate(1).Label, '.TempMatch.mat']);
            TempBout{j} = Temp;
            [Pks, Locs] = findpeaks(conv(Temp.Bout.MaxBoutSeqMatch, SmoothingKernel, 'same'), 'MINPEAKHEIGHT', MPH);
            Pks = Pks(:);
            Locs = Locs(:);
        
            SyllPeaks{j} = [SyllPeaks{j}; [Pks Temp.Bout.T(Locs)' ones(size(Pks))*SongFileNo]];
        end
        if (i == 1)
            for j = 1:length(handles.SyllTemplates.SyllableTemplates),
                [MaxPk, MaxPkTime] = max(conv(TempBout{j}.Bout.MaxBoutSeqMatch, SmoothingKernel, 'same'));
                MaxPkTime = TempBout{j}.Bout.T(MaxPkTime);
                SyllMaxPeaks{j}(SongFileNo,:) = [MaxPk MaxPkTime];
                for k = 1:length(handles.SyllTemplates.SyllableTemplates),
                    if (k == j)
                        SyllOtherPeaks{j}{k}(SongFileNo,:) = SyllMaxPeaks{j}(SongFileNo,:);
                    else
                        MaxPkTimeWindow = find((TempBout{k}.Bout.T(1:length(TempBout{k}.Bout.MaxBoutSeqMatch)) >= (MaxPkTime - 0.010)) & (TempBout{k}.Bout.T(1:length(TempBout{k}.Bout.MaxBoutSeqMatch)) <= (MaxPkTime + 0.010)));
                        TempSmoothedSeqMatch = conv(TempBout{k}.Bout.MaxBoutSeqMatch, SmoothingKernel, 'same');
                        if (~isempty(MaxPkTimeWindow))
                            [MaxOtherPk, MaxOtherPkTime] = max(TempSmoothedSeqMatch(MaxPkTimeWindow));
                            MaxOtherPkTime = TempBout{k}.Bout.T(MaxPkTimeWindow(MaxOtherPkTime));
                            SyllOtherPeaks{j}{k}(SongFileNo,:) = [MaxOtherPk MaxOtherPkTime];
                        else
                            SyllOtherPeaks{j}{k}(SongFileNo,:) = [NaN NaN];
                        end
                    end
                end
            end
        end
    end
    
    Results(i).UnDirSongFiles = SongFiles;
    Results(i).UnDirFileLens = UnDirFileLens;
    
    Results(i).UnDirPeaks{1} = MotifPeaks;
    
    if (exist('SyllPeaks', 'var'))
        for j = 1:length(SyllPeaks),
            Results(i).UnDirPeaks{j+1} = SyllPeaks{j};
        end
        for j = 1:length(SyllPeaks),
            Results(i).UnDirSyllMaxPeaks{j} = SyllMaxPeaks{j};
            Results(i).UnDirSyllMaxOtherPeaks{j} = SyllOtherPeaks{j};
        end
    end
end

% Now analyze the data from the random song comparisons

Fid = fopen('/data2/raghav/DataFromHardDrive/HVC_Microlesions/RandomComparisonSongs/WavSongFiles.txt', 'r');
Temp = textscan(Fid, '%s', 'delimiter', '\n');
fclose(Fid);

SongFiles = Temp{1};

Fid = fopen('/data2/raghav/DataFromHardDrive/HVC_Microlesions/RandomComparisonSongs/ObsSongFiles.txt', 'r');
Temp = textscan(Fid, '%s', 'delimiter', '\n');
fclose(Fid);

SongFiles = [SongFiles; Temp{1}];

Fid = fopen('/data2/raghav/DataFromHardDrive/HVC_Microlesions/RandomComparisonSongs/OKrankSongFiles.txt', 'r');
Temp = textscan(Fid, '%s', 'delimiter', '\n');
fclose(Fid);

SongFiles = [SongFiles; Temp{1}];

MotifPeaks = [];
MotifMaxPeaks = [];
for j = 1:length(handles.SyllTemplates.SyllableTemplates),
    SyllPeaks{j} = [];
    SyllMaxPeaks{j} = [];
end

for SongFileNo = 1:length(SongFiles),
    SongFile = SongFiles{SongFileNo};
    Slash = find((SongFile == '/') | (SongFile == '\'));
    if (~isempty(Slash))
        SongFile = SongFile(Slash(end)+1:end);
    end
    
    if (strfind(SongFile, handles.BirdName))
        continue;
    end
    
    ContinueOrNot = 0;
    if (isfield(handles, 'ExcludeBirds'))
        for i = 1:length(handles.ExcludeBirds),
            if (strfind(SongFile, handles.ExcludeBirds{i}))
                ContinueOrNot = 1;
                break;
            end
        end
    end
    
    if (ContinueOrNot == 1)
        continue;
    end
    
    cd(handles.UnDirRandomSongComparisons_MotifTemplateMatchOutputFilesDir);
    Temp = load([SongFile, '.Motif.TempMatch.mat']);
    [Pks, Locs] = findpeaks(conv(Temp.Bout.MaxBoutSeqMatch, SmoothingKernel, 'same'), 'MINPEAKHEIGHT', MPH);
    Pks = Pks(:);
    Locs = Locs(:);
    [MaxPk, MaxLoc] = max(conv(Temp.Bout.MaxBoutSeqMatch, SmoothingKernel, 'same'));

    MotifPeaks = [MotifPeaks; [Pks Temp.Bout.T(Locs)' ones(size(Pks))*SongFileNo]];
    MotifMaxPeaks = [MotifMaxPeaks; [MaxPk, Temp.Bout.T(MaxLoc) SongFileNo]];
    
    cd(handles.UnDirRandomSongComparisons_SyllTemplateMatchOutputFilesDir);
    for j = 1:length(handles.SyllTemplates.SyllableTemplates),
        Temp = load([SongFile, '.', handles.SyllTemplates.SyllableTemplates{j}{1}.MotifTemplate(1).Label, '.TempMatch.mat']);
        [Pks, Locs] = findpeaks(conv(Temp.Bout.MaxBoutSeqMatch, SmoothingKernel, 'same'), 'MINPEAKHEIGHT', MPH);
        Pks = Pks(:);
        Locs = Locs(:);
        [MaxPk, MaxLoc] = max(conv(Temp.Bout.MaxBoutSeqMatch, SmoothingKernel, 'same'));

        SyllPeaks{j} = [SyllPeaks{j}; [Pks Temp.Bout.T(Locs)' ones(size(Pks))*SongFileNo]];
        SyllMaxPeaks{j} = [SyllMaxPeaks{j}; [MaxPk Temp.Bout.T(MaxLoc) SongFileNo]];
    end
end
    
Results(1).RandomSongComparisonFiles = SongFiles;
Results(1).RandomSongComparisonUnDirPeaks{1} = MotifPeaks;
Results(1).RandomSongComparisonUnDirMaxPeaks{1} = MotifMaxPeaks;

if (exist('SyllPeaks', 'var'))
    for j = 1:length(SyllPeaks),
        Results(1).RandomSongComparisonUnDirPeaks{j+1} = SyllPeaks{j};
        Results(1).RandomSongComparisonUnDirMaxPeaks{j+1} = SyllMaxPeaks{j};
    end
end

% Now analyze the data from the random sound comparisons

Fid = fopen('/home/raghav/RandomSounds/RandomSounds.txt', 'r');
Temp = textscan(Fid, '%s', 'delimiter', '\n');
fclose(Fid);

SongFiles = Temp{1};

MotifPeaks = [];
MotifMaxPeaks = [];
for j = 1:length(handles.SyllTemplates.SyllableTemplates),
    SyllPeaks{j} = [];
    SyllMaxPeaks{j} = [];
end

for SongFileNo = 1:length(SongFiles),
    SongFile = SongFiles{SongFileNo};
    Slash = find((SongFile == '/') | (SongFile == '\'));
    if (~isempty(Slash))
        SongFile = SongFile(Slash(end)+1:end);
    end
    
    cd(handles.UnDirRandomSoundComparison_MotifTemplateMatchOutputFilesDir);
    Temp = load([SongFile, '.Motif.TempMatch.mat']);
    [Pks, Locs] = findpeaks(conv(Temp.Bout.MaxBoutSeqMatch, SmoothingKernel, 'same'), 'MINPEAKHEIGHT', MPH);
    Pks = Pks(:);
    Locs = Locs(:);
    [MaxPk, MaxLoc] = max(conv(Temp.Bout.MaxBoutSeqMatch, SmoothingKernel, 'same'));
    
    MotifPeaks = [MotifPeaks; [Pks Temp.Bout.T(Locs)' ones(size(Pks))*SongFileNo]];
    MotifMaxPeaks = [MotifMaxPeaks; [MaxPk Temp.Bout.T(MaxLoc) SongFileNo]];

    cd(handles.UnDirRandomSoundComparison_SyllTemplateMatchOutputFilesDir);
    for j = 1:length(handles.SyllTemplates.SyllableTemplates),
        Temp = load([SongFile, '.', handles.SyllTemplates.SyllableTemplates{j}{1}.MotifTemplate(1).Label, '.TempMatch.mat']);
        [Pks, Locs] = findpeaks(conv(Temp.Bout.MaxBoutSeqMatch, SmoothingKernel, 'same'), 'MINPEAKHEIGHT', MPH);
        Pks = Pks(:);
        Locs = Locs(:);
        [MaxPk, MaxLoc] = max(conv(Temp.Bout.MaxBoutSeqMatch, SmoothingKernel, 'same'));

        SyllPeaks{j} = [SyllPeaks{j}; [Pks Temp.Bout.T(Locs)' ones(size(Pks))*SongFileNo]];
        SyllMaxPeaks{j} = [SyllMaxPeaks{j}; [MaxPk Temp.Bout.T(MaxLoc) SongFileNo]];
    end
end
    
Results(1).RandomSoundComparisonFiles = SongFiles;
Results(1).RandomSoundComparisonUnDirPeaks{1} = MotifPeaks;
Results(1).RandomSoundComparisonUnDirMaxPeaks{1} = MotifMaxPeaks;

if (exist('SyllPeaks', 'var'))
    for j = 1:length(SyllPeaks),
        Results(1).RandomSoundComparisonUnDirPeaks{j+1} = SyllPeaks{j};
        Results(1).RandomSoundComparisonUnDirMaxPeaks{j+1} = SyllMaxPeaks{j};
    end
end

% Now analyze the shuffled song comparisons
Fid = fopen(handles.UnDirSongFileList{1}, 'r');
Temp = textscan(Fid, '%s', 'delimiter', '\n');
fclose(Fid);

SongFiles = Temp{1}(1:25);

MotifPeaks = [];
for j = 1:length(handles.SyllTemplates.SyllableTemplates),
    SyllPeaks{j} = [];
end

for SongFileNo = 1:length(SongFiles),
    SongFile = SongFiles{SongFileNo};
    Slash = find((SongFile == '/') | (SongFile == '\'));
    if (~isempty(Slash))
        SongFile = SongFile(Slash(end)+1:end);
    end
    cd(handles.UnDirShuffledSongComparisons_MotifTemplateMatchOutputFilesDir);
    for Repetition = 1:100,
        Temp = load([SongFile, '.Motif.', num2str(Repetition), '.TempMatch.mat']);
        [Pks, Locs] = findpeaks(conv(Temp.Bout.MaxBoutSeqMatch, SmoothingKernel, 'same'), 'MINPEAKHEIGHT', MPH);
        Pks = Pks(:);
        Locs = Locs(:);
        [MaxPk, MaxLoc] = max(conv(Temp.Bout.MaxBoutSeqMatch, SmoothingKernel, 'same'));
    
        MotifPeaks = [MotifPeaks; [Pks Temp.Bout.T(Locs)' ones(size(Pks))*SongFileNo]];
        MotifMaxPeaks = [MotifMaxPeaks; [MaxPk Temp.Bout.T(MaxLoc) SongFileNo]];
    end
    
    cd(handles.UnDirShuffledSongComparisons_SyllTemplateMatchOutputFilesDir);
    for j = 1:length(handles.SyllTemplates.SyllableTemplates),
        for Repetition = 1:100,
            Temp = load([SongFile, '.', handles.SyllTemplates.SyllableTemplates{j}{1}.MotifTemplate(1).Label, '.', num2str(Repetition), '.TempMatch.mat']);
            [Pks, Locs] = findpeaks(conv(Temp.Bout.MaxBoutSeqMatch, SmoothingKernel, 'same'), 'MINPEAKHEIGHT', MPH);
            Pks = Pks(:);
            Locs = Locs(:);
            [MaxPk, MaxLoc] = max(conv(Temp.Bout.MaxBoutSeqMatch, SmoothingKernel, 'same'));

            SyllPeaks{j} = [SyllPeaks{j}; [Pks Temp.Bout.T(Locs)' ones(size(Pks))*SongFileNo]];
            SyllMaxPeaks{j} = [SyllMaxPeaks{j}; [MaxPk Temp.Bout.T(MaxLoc) SongFileNo]];
        end
    end
end

Results(1).ShuffledSongComparisonsSongFiles = SongFiles;
Results(1).ShuffledSongComparisonsUnDirPeaks{1} = MotifPeaks;
Results(1).ShuffledSongComparisonsUnDirMaxPeaks{1} = MotifMaxPeaks;

if (exist('SyllPeaks', 'var'))
    for j = 1:length(SyllPeaks),
        Results(1).ShuffledSongComparisonsUnDirPeaks{j+1} = SyllPeaks{j};
        Results(1).ShuffledSongComparisonsUnDirMaxPeaks{j+1} = SyllMaxPeaks{j};
    end
end

disp('Finished analysis');