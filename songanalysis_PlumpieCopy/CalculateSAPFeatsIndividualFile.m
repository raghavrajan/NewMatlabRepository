function [Feats, EntireAmplitude] = CalculateSAPFeatsIndividualFile(DirectoryName, NoteDir, FileName, FileType)

PresentDir = pwd;
SongFile = FileName;
SlashIndex = find(SongFile == '/');
if (~isempty(SlashIndex))
    SongFile = SongFile(SlashIndex(end) + 1:end);
else
    SongFile = SongFile;
end

SyllNo = 1;
disp(SongFile);
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
            end
        end
    end
    Song = resample(Song, 44100, Fs);
    Fs = 44100;
    Time = (1:1:length(Song))/Fs;
    NoteFile = load([NoteDir, '/', SongFile, '.not.mat']);
    for i = 1:length(NoteFile.labels),
        SyllableStart = find(Time >= NoteFile.onsets(i)/1000, 1, 'first');
        if ((SyllableStart - round(Fs * 0.5)) <= 0)
            SyllableStart = 1;
        else
            SyllableStart = SyllableStart - round(Fs * 0.5);
        end

        SyllableEnd = find(Time > NoteFile.offsets(i)/1000, 1, 'first');
        if ((SyllableEnd + round(Fs * 0.5)) > length(Song))
            SyllableEnd = length(Song);
        else
            SyllableEnd = SyllableEnd + round(Fs * 0.5);
        end

        [m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude , m_Freq, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight] = deriv(Song(SyllableStart:SyllableEnd), Fs);

        T = linspace(Time(SyllableStart), Time(SyllableEnd), length(m_Entropy));        
        StartIndex = find((T>=NoteFile.onsets(i)/1000), 1, 'first');
        EndIndex = find((T>NoteFile.offsets(i)/1000), 1, 'first');
        Feats.AM(SyllNo) = mean(m_AM(StartIndex:EndIndex));
        Feats.Entropy(SyllNo) = mean(m_Entropy(StartIndex:EndIndex));
        Feats.Freq(SyllNo) = mean(m_Freq(StartIndex:EndIndex));
        Feats.Pitch(SyllNo) = mean(Pitch_chose(StartIndex:EndIndex));
        Feats.FM(SyllNo) = mean(m_FM(StartIndex:EndIndex));
        Feats.Amp(SyllNo) = mean(m_amplitude(StartIndex:EndIndex));
        Feats.Freq(SyllNo) = mean(m_Freq(StartIndex:EndIndex));
        Feats.PG(SyllNo) = mean(m_PitchGoodness(StartIndex:EndIndex));
        Feats.Dur(SyllNo) = NoteFile.offsets(i) - NoteFile.onsets(i);
        EntireAmplitude{SyllNo} = m_amplitude(StartIndex:EndIndex);
        SyllNo = SyllNo + 1;
    end
catch
    disp('Skipped file');
end
