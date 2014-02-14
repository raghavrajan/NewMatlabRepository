function [Notes, ANotes] = INALoadNoteFiles(DataDir, RecFileDir, NoteFiles, NoteFileDir, FileType, ContinuousOrNot)

PresentDir = pwd;
cd(NoteFileDir);

if (strfind(ContinuousOrNot, 'continuous'))
    CNotes.Onsets{1} = [];
    CNotes.Offsets{1} = [];
    CNotes.Labels{1} = [];
else
    Notes = [];
end


FTime = 0;
for i = 1:length(NoteFiles),
    Temp = load(NoteFiles(i).name);
    disp(['Loaded ', NoteFiles(i).name]);
    
    Index = strfind(NoteFiles(i).name, '.not.mat');
    SongFile = NoteFiles(i).name(1:Index-1);
    
    if (strfind(FileType, 'okrank'))
        [Song, Fs] = SSAReadOKrankData(DataDir, RecFileDir, SongFile, 1);
    else
        if (strfind(FileType, 'wav'))
            cd(DataDir);
            [Song, Fs] = wavread(SongFile);
            cd(NoteFileDir);
        else
            if (strfind(FileType, 'obs'));
                [Song, Fs] = SSASoundIn(DataDir, RecFileDir, SongFile, 'obs0r');
                Song = Song * 5/32768;
            end
        end
    end
    
    Time = (0:1:(length(Song)-1))/Fs;
    
    if (strfind(ContinuousOrNot, 'continuous'))
        CNotes.Onsets{1} = [CNotes.Onsets{1}; (Temp.onsets/1000 + FTime)];
        CNotes.Offsets{1} = [CNotes.Offsets{1}; (Temp.offsets/1000 + FTime)];
        CNotes.Labels{1} = [CNotes.Labels{1} (Temp.labels)];
        FTime = FTime + Time(end);
    end
    Notes.Onsets{i} = Temp.onsets/1000;
    Notes.Offsets{i} = Temp.offsets/1000;
    Notes.Labels{i} = Temp.labels;
    Notes.FileLength{i} = Time(end);
    Notes.FileName{i} = SongFile;
    Notes.FileType{i} = FileType;
    Notes.DataDir{i} = DataDir;
    Notes.RecFileDir{i} = RecFileDir;
end

if (strfind(ContinuousOrNot, 'continuous'))
    Indices = find((Notes.Labels{1} >= 65) & (Notes.Labels{1} <= 90));
    if (~isempty(Indices))
        i = 1;
        while (i <= length(Indices)),
            if (i ~= length(Indices))
                if (Indices(i+1) == (Indices(i) + 1))
                    Notes.Labels{1}(Indices(i)) = lower(Notes.Labels{1}(Indices(i)));
                    Notes.Offsets{1}(Indices(i)) = Notes.Offsets{1}(Indices(i+1));
                    Notes.Onsets{1}(Indices(i+1)) = -1000;
                    i = i + 1;
                else
                    Notes.Labels{1}(Indices(i)) = lower(Notes.Labels{1}(Indices(i)));
                end
            end
            i = i + 1;
        end
    end
    Indices = find(Notes.Onsets{1} < 0);
    Notes.Onsets{1}(Indices) = [];
    Notes.Offsets{1}(Indices) = [];
    Notes.Labels{1}(Indices) = [];
end

ANotes = Notes;

if (strfind(ContinuousOrNot, 'continuous'))
    Notes = CNotes;
end

disp('Finished loading all note files');
cd(PresentDir);
