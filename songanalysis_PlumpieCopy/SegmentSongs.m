function [NoteInfo] = SegmentSongs(DirectoryName, FileList, FileType, Threshold, PlotOption)

min_int = 0.002;
min_dur = 0.015;
sm_win = 0.008;

FFTWinSize = sm_win; % in sec
FFTWinOverlap = 0.95;

if (DirectoryName(end) ~= '/')
    DirectoryName(end + 1) = '/';
end

cd(DirectoryName);

NoteInfo.FileName = [];
NoteInfo.Onsets = [];
NoteInfo.Offsets = [];
NoteInfo.AM = [];
NoteInfo.FM = [];
NoteInfo.Entropy = [];
NoteInfo.Amp = [];
NoteInfo.Freq = [];
NoteInfo.PG = [];
NoteInfo.Pitch = [];
NoteInfo.Duration = [];

Fid = fopen(FileList, 'r');
SongFile = fgetl(Fid);

while (ischar(SongFile(1)))
    Slash = find(SongFile == '/');
    if (~isempty(Slash))
        SongFile = SongFile(Slash(end)+1:end);
    end
    disp(SongFile);
    try
        if (strfind(FileType, 'okrank'))
            [Song, Fs] = ReadOKrankData(DirectoryName, SongFile, 1);
        else
            if (strfind(FileType, 'wav'))
                [Song, Fs] = wavread(SongFile);
            else
                if (strfind(FileType, 'obs'))
                    channel_string = strcat('obs',num2str(0),'r');
                    [Song, Fs] = soundin_copy(DirectoryName, SongFile, channel_string);
    
                    % Convert to uV - 5V on the data acquisition is 32768
                    Song = Song * 5/32768;
                end
            end
        end
    catch
        continue;
    end
    if (isempty(Song))
        continue;
    end
    [m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude , m_Freq, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight ]=deriv(Song, Fs);
    
    WinSize = round(FFTWinSize * Fs);
    WinOverlap = round(FFTWinOverlap * WinSize);

    [S, F, T, P] = spectrogram(Song, hamming(WinSize), WinOverlap, WinSize, Fs);
    
    Freq = find((F >= 2000) & (F <= 7000));
    Power = sum(log10(abs(S(Freq,:))));
    % Power = 10*log10(sum(P));
    ThreshCrossing = Power > Threshold;
    [onsets, offsets] = segment_song(Power, (1/(T(2) - T(1))), min_int*1000, min_dur*1000, Threshold);
    onsets = onsets/1000;
    offsets = offsets/1000;
%     h = [1 -1];
%     Temp = zeros(size(ThreshCrossing));
%     Temp(ThreshCrossing > 0) = 1;
%     Transitions = conv(Temp, h, 'same');
%     onsets = T(Transitions > 0);
%     offsets = T(Transitions < 0);
%     if (length(offsets) > length(onsets))
%         offsets(1) = [];
%     else
%         if (length(offsets) < length(onsets))
%             onsets(end) = [];
%         end
%     end
%     
%     ShortNotes = find((offsets - onsets) < 0.015);
%     onsets(ShortNotes) = [];
%     offsets(ShortNotes) = [];
    
    if (size(onsets,1) < size(onsets, 2))
        onsets = onsets';
    end
    if (size(offsets,1) < size(offsets, 2))
        offsets = offsets';
    end
    
    if (strfind(PlotOption, 'on'))
        PlotSpectrogram(DirectoryName, SongFile, FileType);
        hold on;
        plot(T, (Power + 70)*100);
        for j = 1:length(onsets),
            plot([onsets(j) onsets(j) offsets(j) offsets(j)], [0 8000 8000 0],'r');
        end
        zoom xon;
        uiwait(gcf);
    end
    Time = linspace(T(1), T(end), length(m_AM));
    for NoteIndex = 1:length(onsets),
        StartI = find(Time < onsets(NoteIndex), 1, 'last');
        StartI = StartI + 1;
        EndI = find(Time > offsets(NoteIndex), 1, 'first');
        EndI = EndI - 1;
        FileName{NoteIndex} = SongFile;
        AM(NoteIndex,:) = mean(m_AM(StartI:EndI));
        FM(NoteIndex,:) = mean(m_FM(StartI:EndI));
        Entropy(NoteIndex,:) = mean(m_Entropy(StartI:EndI));
        Amp(NoteIndex,:) = mean(m_amplitude(StartI:EndI));
        Freq(NoteIndex,:) = mean(m_Freq(StartI:EndI));
        PG(NoteIndex,:) = mean(m_PitchGoodness(StartI:EndI));
        Pitch(NoteIndex,:) = mean(Pitch_chose(StartI:EndI));
        Duration(NoteIndex,:) = offsets(NoteIndex) - onsets(NoteIndex);
    end
    
    NoteInfo.FileName = [NoteInfo.FileName FileName];
    NoteInfo.Onsets = [NoteInfo.Onsets; onsets];
    NoteInfo.Offsets = [NoteInfo.Offsets; offsets];
    NoteInfo.AM = [NoteInfo.AM; AM];
    NoteInfo.FM = [NoteInfo.FM; FM];
    NoteInfo.Entropy = [NoteInfo.Entropy; Entropy];
    NoteInfo.Amp = [NoteInfo.Amp; Amp];
    NoteInfo.Freq = [NoteInfo.Freq; Freq];
    NoteInfo.PG = [NoteInfo.PG; PG];
    NoteInfo.Pitch = [NoteInfo.Pitch; Pitch];
    NoteInfo.Duration = [NoteInfo.Duration; Duration];
    clear onsets offsets AM FM Entropy Amp Freq PG Pitch Duration FileName;
    SongFile = fgetl(Fid);
end
fclose(Fid);
