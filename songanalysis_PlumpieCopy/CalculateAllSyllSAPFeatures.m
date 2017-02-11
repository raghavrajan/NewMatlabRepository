function [SyllFeats, MotifLengths] = CalculateAllSyllSAPFeatures(NoteFileList, RawDataDirectory, FileType, Motif)

PresentDir = pwd;
Fid = fopen(NoteFileList, 'r');
Files = textscan(Fid, '%s', 'DeLimiter', '\n');
Files = Files{1};

if (ispc)
    if (RawDataDirectory(end) ~= '\')
        RawDataDirectory(end+1) = '\';
    end
else
    if (RawDataDirectory(end) ~= '/')
        RawDataDirectory(end+1) = '/';
    end
end

AllFeats = [];
AllLabels = [];
MotifLengths = [];

for i = 1:length(Files),
    if (~exist([Files{i}, '.not.mat'], 'file'))
        break;
    end
    Notes{i} = load([Files{i}, '.not.mat']);
    MotifMatches = strfind(Notes{i}.labels, Motif);
    TempMotifLengths = Notes{i}.offsets(MotifMatches + length(Motif) - 1)/1000 - Notes{i}.onsets(MotifMatches)/1000;
    TempMotifLengths = TempMotifLengths(:);
    MotifLengths = [MotifLengths; TempMotifLengths];
    
    if (strfind(FileType,'obs'))
        [Song,Fs] = soundin_copy(RawDataDirectory, Files{i}, 'obs0r');
        Song = Song * 1/32768;
    else
        if (strfind(FileType,'wav'));
            cd(RawDataDirectory);
            disp(Files{i});
            [Song, Fs] = wavread(Files{i});
            cd(PresentDir);
        else 
            if (strfind(FileType, 'okrank'))
                [Song, Fs] = ReadOKrankData(RawDataDirectory, Files{i}, 1);
                Song = Song/10;
            end
        end
    end
    Time = [1:1:length(Song)]/Fs;

    [TempFeats] = CalculateSAPFeatsWithOnsets(Song, Time, Fs, Notes{i}.onsets/1000, Notes{i}.offsets/1000);    
    AllFeats = [AllFeats; TempFeats];
    AllLabels = [AllLabels Notes{i}.labels];
end

UniqueLabels = unique(AllLabels);

for i = 1:length(UniqueLabels),
    Matches = find(AllLabels == UniqueLabels(i));
    SyllFeats{i}.Label = UniqueLabels(i);
    SyllFeats{i}.Feats = AllFeats(Matches,:);
end

disp('Finished');