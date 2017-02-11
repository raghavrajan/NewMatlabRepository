function [Feats] = CalculateSAPFeatsAllSyllables(DirectoryName, NoteDir, FileList, FileType, MaxFiles)

PresentDir = pwd;
Fid = fopen(FileList, 'r');
SongFile = fgetl(Fid);
SlashIndex = find(SongFile == '/');
if (~isempty(SlashIndex))
    SongFile = SongFile(SlashIndex(end) + 1:end);
else
    SongFile = SongFile;
end

SyllNo = 1;
FileIndex = 0;
while (ischar(SongFile(1)))
    FileIndex = FileIndex + 1;
    FileSyllNo = 1;
    
    if (FileIndex > MaxFiles)
        break;
    end
    
    SlashIndex = find((SongFile == '/') | (SongFile == '\'));
    if (~isempty(SlashIndex))
        SongFile = SongFile(SlashIndex(end)+1:end);
    end
    disp(SongFile);
    
    Feats.FileData{FileIndex}.FileName = SongFile;
    
    try
        if (strfind(FileType, 'okrank'))
            [Song, Fs] = ReadOKrankData(DirectoryName, SongFile, 1);
            Song = Song/10;
        else
            if (strfind(FileType, 'obs'))
                if ispc
                    [Song, Fs] = soundin_copy([DirectoryName, '\'], SongFile, 'obs0r');
                else
                    [Song, Fs] = soundin_copy([DirectoryName, '/'], SongFile, 'obs0r');
                end
                Song = Song/32768;
            else
                if (strfind(FileType, 'wav'))
                    cd(DirectoryName);
                    [Song, Fs] = wavread(SongFile);
                    cd(PresentDir);
                else
                    if (strfind(FileType, 'filt'))
                        [Song, Fs] = soundin_copy(DirectoryName, [SongFile, '.filt'], 'filt');
                    end
                end
            end
        end
        Song = resample(Song, 44100, Fs);
        Fs = 44100;
        Time = (1:1:length(Song))/Fs;
        NoteFile = load([NoteDir, '/', SongFile, '.not.mat']);
        Sylls = 1:1:length(NoteFile.labels);
        for i = 1:length(Sylls),
            SyllableStart = find(Time >= NoteFile.onsets(Sylls(i))/1000, 1, 'first');
            if ((SyllableStart - round(Fs * 0.5)) <= 0)
                SyllableStart = 1;
            else
                SyllableStart = SyllableStart - round(Fs * 0.5);
            end

            SyllableEnd = find(Time > NoteFile.offsets(Sylls(i))/1000, 1, 'first');
            if ((SyllableEnd + round(Fs * 0.5)) > length(Song))
                SyllableEnd = length(Song);
            else
                SyllableEnd = SyllableEnd + round(Fs * 0.5);
            end

            [m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude , m_Freq, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight] = deriv(Song(SyllableStart:SyllableEnd), Fs);

            T = linspace(Time(SyllableStart), Time(SyllableEnd), length(m_Entropy));        
            
           
            % This is the start and end index for the syllable to pull out
            % the SAP values associated with the syllable
            
            StartIndex = find((T>=NoteFile.onsets(Sylls(i))/1000), 1, 'first');
            EndIndex = find((T>NoteFile.offsets(Sylls(i))/1000), 1, 'first');
            
            
            Feats.AM(SyllNo) = mean(m_AM(StartIndex:EndIndex));
            Feats.Entropy(SyllNo) = mean(m_Entropy(StartIndex:EndIndex));
            % Feats.Pitch(SyllNo) = mean(Pitch_chose(StartIndex:EndIndex)); % harmonic pitch
            Feats.FM(SyllNo) = mean(m_FM(StartIndex:EndIndex));
            Feats.Amp(SyllNo) = mean(m_amplitude(StartIndex:EndIndex));
            Feats.Freq(SyllNo) = mean(m_Freq(StartIndex:EndIndex)); % mean frequency
            Feats.PG(SyllNo) = mean(m_PitchGoodness(StartIndex:EndIndex));
            Feats.Dur(SyllNo) = NoteFile.offsets(Sylls(i)) - NoteFile.onsets(Sylls(i)); % in ms
            Feats.EntropyVar(SyllNo) = var(m_Entropy(StartIndex:EndIndex));
            Feats.Label(SyllNo) = NoteFile.labels(Sylls(i));
            
            Feats.FileData{FileIndex}.AM(FileSyllNo) = mean(m_AM(StartIndex:EndIndex));
            Feats.FileData{FileIndex}.Entropy(FileSyllNo) = mean(m_Entropy(StartIndex:EndIndex));
            Feats.FileData{FileIndex}.Freq(FileSyllNo) = mean(m_Freq(StartIndex:EndIndex));
            Feats.FileData{FileIndex}.Pitch(FileSyllNo) = mean(Pitch_chose(StartIndex:EndIndex));
            Feats.FileData{FileIndex}.FM(FileSyllNo) = mean(m_FM(StartIndex:EndIndex));
            Feats.FileData{FileIndex}.Amp(FileSyllNo) = mean(m_amplitude(StartIndex:EndIndex));
            Feats.FileData{FileIndex}.Freq(FileSyllNo) = mean(m_Freq(StartIndex:EndIndex));
            Feats.FileData{FileIndex}.PG(FileSyllNo) = mean(m_PitchGoodness(StartIndex:EndIndex));
            Feats.FileData{FileIndex}.Dur(FileSyllNo) = NoteFile.offsets(Sylls(i)) - NoteFile.onsets(Sylls(i));
            Feats.FileData{FileIndex}.EntropyVar(FileSyllNo) = var(m_Entropy(StartIndex:EndIndex));
            Feats.FileData{FileIndex}.Label(FileSyllNo) = NoteFile.labels(Sylls(i));
            FileSyllNo = FileSyllNo + 1;
            SyllNo = SyllNo + 1;
        end
    catch
        disp('Skipped file');
    end
    SongFile = fgetl(Fid);
    SlashIndex = find(SongFile == '/');
    if (~isempty(SlashIndex))
        SongFile = SongFile(SlashIndex(end) + 1:end);
    else
        SongFile = SongFile;
    end
end        
fclose(Fid);