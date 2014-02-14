function [BoutStats] = CalcBoutIntroNoteDist(RawDataDir, FileType, NoteDirName, NoteFiles, Motif, Motif2, IntroSyllables, ContinuousOrNot)

InterBoutInterval = 0.5; % in seconds

if ispc
    if (RawDataDir(end) ~= '\')
        RawDataDir(end+1) = '\';
    end
else
    if (RawDataDir(end) ~= '/')
        RawDataDir(end+1) = '/';
    end
end
      
BoutNo = 1;

cd(NoteDirName);

PrevFTime = 0;
LastSongBoutTime = -100;


if (strfind(ContinuousOrNot, 'continuous'))
    Labels = [];
    Onsets = [];
    Offsets = [];
    for NoteFileNo = 1:length(NoteFiles),
        Notes = load(NoteFiles(NoteFileNo).name);
        disp(NoteFiles(NoteFileNo).name);
        disp(Notes.labels);
        Notes.onsets = Notes.onsets/1000;
        Notes.offsets = Notes.offsets/1000;
        
        Index = strfind(NoteFiles(NoteFileNo).name, '.not.mat');
        SongFile = NoteFiles(NoteFileNo).name(1:Index-1);

        if (strfind(FileType, 'okrank'))
            [Song, Fs] = ReadOKrankData(RawDataDir, SongFile, 1);
        else
            if (strfind(FileType, 'wav'))
                cd(RawDataDir);
                [Song, Fs] = wavread(SongFile);
                cd(NoteDirName);
            else
                if (strfind(FileType, 'obs'));
                    [Song, Fs] = soundin_copy(RawDataDir, SongFile, 'obs0r');
                    Song = Song * 5/32768;
                end
            end
        end

        Time = (0:1:(length(Song)-1))/Fs;

        Labels = [Labels Notes.labels];
        Onsets = [Onsets; (Notes.onsets + PrevFTime)];
        Offsets = [Offsets; (Notes.offsets + PrevFTime)];
        
        FTime = PrevFTime + length(Song)/Fs;
        PrevFTime = FTime;
    end
    Indices = find((Labels >= 65) & (Labels <= 90));
    if (~isempty(Indices))
        i = 1;
        while (i <= length(Indices)),
            if (i ~= length(Indices))
                if (Indices(i+1) == (Indices(i) + 1))
                    Labels(Indices(i)) = lower(Labels(Indices(i)));
                    Offsets(Indices(i)) = Offsets(Indices(i+1));
                    Onsets(Indices(i+1)) = -1000;
                    i = i + 1;
                else
                    Labels(Indices(i)) = lower(Labels(Indices(i)));
                end
            end
            i = i + 1;
        end
    end
    Indices = find(Onsets < 0);
    Onsets(Indices) = [];
    Offsets(Indices) = [];
    Labels(Indices) = [];
    
    Bouts = find((Onsets(2:end) - Offsets(1:(end-1))) > InterBoutInterval);

    for i = 0:length(Bouts),
        if (i == 0)
            StartIndex = 1;
        else
            StartIndex = Bouts(i) + 1;
        end

        if (i == length(Bouts))
            EndIndex = length(Onsets);
        else
            EndIndex = Bouts(i+1);
        end

        Motifs = union(strfind(Labels(StartIndex:EndIndex), Motif), strfind(Labels(StartIndex:EndIndex), Motif2));
        if (~isempty(Motifs))
            Motifs = Motifs + StartIndex - 1;
            if ((Onsets(StartIndex) > InterBoutInterval) && (Offsets(EndIndex) < (PrevFTime - InterBoutInterval)))
                BoutStats.IntroNoteNos(BoutNo) = 0;
                for j = 1:length(IntroSyllables),
                    BoutStats.IntroNoteNos(BoutNo) = BoutStats.IntroNoteNos(BoutNo) + length(find(Labels(StartIndex:Motifs(1)) == IntroSyllables(j)));
                end
                BoutStats.BoutLength(BoutNo) = Offsets(EndIndex) - Onsets(StartIndex);
                if (StartIndex == 1)
                    BoutStats.VocalInterval(BoutNo) = Onsets(1);
                else
                    BoutStats.VocalInterval(BoutNo) = Onsets(StartIndex) - Offsets(StartIndex - 1);
                end
                BoutStats.Onset(BoutNo) = Onsets(StartIndex);
                BoutStats.BoutInterval(BoutNo) = Onsets(StartIndex) - LastSongBoutTime;
                LastSongBoutTime = Offsets(EndIndex);
                BoutNo = BoutNo + 1;
            end
        end
    end
else
    for NoteFileNo = 1:length(NoteFiles),
        Notes = load(NoteFiles(NoteFileNo).name);
        disp(NoteFiles(NoteFileNo).name);
        disp(Notes.labels);
        Notes.onsets = Notes.onsets/1000;
        Notes.offsets = Notes.offsets/1000;

        Index = strfind(NoteFiles(NoteFileNo).name, '.not.mat');
        SongFile = NoteFiles(NoteFileNo).name(1:Index-1);

        if (strfind(FileType, 'okrank'))
            [Song, Fs] = ReadOKrankData(RawDataDir, SongFile, 1);
        else
            if (strfind(FileType, 'wav'))
                cd(RawDataDir);
                [Song, Fs] = wavread(SongFile);
                cd(NoteDirName);
            else
                if (strfind(FileType, 'obs'));
                    [Song, Fs] = soundin_copy(RawDataDir, SongFile, 'obs0r');
                    Song = Song * 5/32768;
                end
            end
        end

        Time = (0:1:(length(Song)-1))/Fs;

        if (isempty(Notes.offsets))
            continue;
        end

        Bouts = find((Notes.onsets(2:end) - Notes.offsets(1:(end-1))) > InterBoutInterval);

        for i = 0:length(Bouts),
            if (i == 0)
                StartIndex = 1;
            else
                StartIndex = Bouts(i) + 1;
            end

            if (i == length(Bouts))
                EndIndex = length(Notes.onsets);
            else
                EndIndex = Bouts(i+1);
            end

            Motifs = union(strfind(Notes.labels(StartIndex:EndIndex), Motif), strfind(Notes.labels(StartIndex:EndIndex), Motif2));
            if (~isempty(Motifs))
                Motifs = Motifs + StartIndex - 1;
                if ((Notes.onsets(StartIndex) > InterBoutInterval) && (Notes.offsets(EndIndex) < (Time(end) - InterBoutInterval)))
                    BoutStats.IntroNoteNos(BoutNo) = 0;
                    for j = 1:length(IntroSyllables),
                        BoutStats.IntroNoteNos(BoutNo) = BoutStats.IntroNoteNos(BoutNo) + length(find(Notes.labels(StartIndex:Motifs(1)) == IntroSyllables(j)));
                    end
                    BoutStats.BoutLength(BoutNo) = Notes.offsets(EndIndex) - Notes.onsets(StartIndex);
                    if (StartIndex == 1)
                        BoutStats.BoutInterval(BoutNo) = Notes.onsets(1) - Time(1);
                    else
                        BoutStats.BoutInterval(BoutNo) = Notes.onsets(StartIndex) - Notes.offsets(StartIndex - 1);
                    end
                    BoutStats.FileName{BoutNo} = SongFile;
                    BoutNo = BoutNo + 1;
                end
            end
        end
    end    
end

disp('Finished analysing notes');