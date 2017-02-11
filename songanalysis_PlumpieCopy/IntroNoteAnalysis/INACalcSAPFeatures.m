function [MotifIntroNoteIndices, NonMotifIntroNoteIndices] = INACalcSAPFeatures(MotifIntroNoteIndices, NonMotifIntroNoteIndices, Onsets, Offsets, SongFile, FileType, DataDir, RecFileDir)

PresentDir = pwd;

if (strfind(FileType, 'okrank'))
    [Song, Fs] = SSAReadOKrankData(DataDir, RecFileDir, SongFile, 1);
else
    disp(SongFile);
    if (strfind(FileType, 'wav'))
        cd(DataDir);
        [Song, Fs] = wavread(SongFile);
        cd(PresentDir);
    else
        if (strfind(FileType, 'obs'));
            [Song, Fs] = SSASoundIn(DataDir, RecFileDir, SongFile, 'obs0r');
            Song = Song * 5/32768;
        end
    end
end
    
Song = resample(Song, 44100, Fs);
Fs = 44100;
Time = (0:1:(length(Song)-1))/Fs;


if (~isfield(MotifIntroNoteIndices, 'Onsets'))
    return;
end
for i = 1:length(MotifIntroNoteIndices.Onsets),
    for j = 1:length(MotifIntroNoteIndices.Onsets{i}),
        SyllableStart = find(Time <= MotifIntroNoteIndices.Onsets{i}(j), 1, 'last');
        if ((SyllableStart - round(Fs * 0.5)) <= 0)
            SyllableStart = 1;
        else
            SyllableStart = SyllableStart - round(Fs * 0.5);
        end
        
        SyllableEnd = find(Time >= MotifIntroNoteIndices.Offsets{i}(j), 1, 'first');
        if ((SyllableEnd + round(Fs * 0.5)) > length(Song))
            SyllableEnd = length(Song);
        else
            SyllableEnd = SyllableEnd + round(Fs * 0.5);
        end
        
        [m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude , m_Freq, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight] = deriv(Song(SyllableStart:SyllableEnd), Fs);
        
        T = linspace(Time(SyllableStart), Time(SyllableEnd), length(m_Entropy));        
        
        StartIndex = find(T <= MotifIntroNoteIndices.Onsets{i}(j), 1, 'last');
        EndIndex = find(T >= MotifIntroNoteIndices.Offsets{i}(j), 1, 'first');
        MotifIntroNoteIndices.PG{i}(j) = mean(m_PitchGoodness(StartIndex:EndIndex));
        MotifIntroNoteIndices.Amplitude{i}(j) = mean(m_amplitude(StartIndex:EndIndex));
        MotifIntroNoteIndices.FM{i}(j) = mean(m_FM(StartIndex:EndIndex));
        MotifIntroNoteIndices.AM{i}(j) = mean(m_AM(StartIndex:EndIndex));
        MotifIntroNoteIndices.Pitch{i}(j) = mean(Pitch_chose(StartIndex:EndIndex));
        MotifIntroNoteIndices.Entropy{i}(j) = mean(m_Entropy(StartIndex:EndIndex));
        MotifIntroNoteIndices.MeanFreq{i}(j) = mean(m_Freq(StartIndex:EndIndex));
    end
end

for i = 1:length(MotifIntroNoteIndices.Onsets),
    if (~isempty(MotifIntroNoteIndices.Onsets{i}))
        MotifIntroNoteIndices.Durs{i} = MotifIntroNoteIndices.Offsets{i} - MotifIntroNoteIndices.Onsets{i};
%        [Val, Ind] = min(abs(onsets - Onsets(MotifIntroNoteIndices.Indices{i}(end) + 1)));
        MotifIntroNoteIndices.Intervals{i} = [MotifIntroNoteIndices.Onsets{i}(2:end); Onsets(MotifIntroNoteIndices.Indices{i}(end) + 1)] - MotifIntroNoteIndices.Offsets{i};
    end
end

for i = 1:length(NonMotifIntroNoteIndices.Onsets),
    for j = 1:length(NonMotifIntroNoteIndices.Onsets{i}),
 %       [Val, Ind] = min(abs(onsets - NonMotifIntroNoteIndices.Onsets{i}(j)));
 %       NonMotifIntroNoteIndices.Onsets{i}(j) = onsets(Ind);
 %       [Val, Ind] = min(abs(offsets - NonMotifIntroNoteIndices.Offsets{i}(j)));
 %       NonMotifIntroNoteIndices.Offsets{i}(j) = offsets(Ind);

        StartIndex = find(T <= NonMotifIntroNoteIndices.Onsets{i}(j), 1, 'last');
        EndIndex = find(T >= NonMotifIntroNoteIndices.Offsets{i}(j), 1, 'first');
        NonMotifIntroNoteIndices.PG{i}(j) = mean(m_PitchGoodness(StartIndex:EndIndex));
        NonMotifIntroNoteIndices.Amplitude{i}(j) = mean(m_amplitude(StartIndex:EndIndex));
        NonMotifIntroNoteIndices.FM{i}(j) = mean(m_FM(StartIndex:EndIndex));
        NonMotifIntroNoteIndices.AM{i}(j) = mean(m_AM(StartIndex:EndIndex));
        NonMotifIntroNoteIndices.Pitch{i}(j) = mean(Pitch_chose(StartIndex:EndIndex));
        NonMotifIntroNoteIndices.Entropy{i}(j) = mean(m_Entropy(StartIndex:EndIndex));
        NonMotifIntroNoteIndices.MeanFreq{i}(j) = mean(m_Freq(StartIndex:EndIndex));
    end
end

for i = 1:length(NonMotifIntroNoteIndices.Onsets),
    if (~isempty(NonMotifIntroNoteIndices.Onsets{i}))
        NonMotifIntroNoteIndices.Durs{i} = NonMotifIntroNoteIndices.Offsets{i} - NonMotifIntroNoteIndices.Onsets{i};
    end
end