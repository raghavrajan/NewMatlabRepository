function [NoteInfo] = DoAutoAnalysis(DirectoryName, FileList, FileType)

cd(DirectoryName);

NoteInfo = [];
for i = 1:length(FileList),
    SongFileName = FileList(i).name;
    if (strfind(FileType, 'okrank'))
        [RawData, Fs] = ReadOKrankData(DirectoryName, SongFileName, 0);
    else
        if (strfind(FileType, 'obs'))
            [RawData, Fs] = soundin_copy(DirectoryName, SongFileName, 'obs0r');   
            RawData = RawData/32768;
        else
            if (strfind(FileType, 'wav'))
                [RawData, Fs] = wavread(SongFileName);    
            end
        end
    end
    Time = (1:1:length(RawData))/Fs;

    % Now using an 8 pole butterworth bandpass filter as default.
    [b,a]=butter(8,[300*2/Fs, 10000*2/Fs]);

    FiltSong=filtfilt(b, a, RawData);
  
    if length(RawData) ~= length(FiltSong) 
        disp(['warning! bandpass: input and output file lengths do not match!']);
    end

    RawData = FiltSong;
    clear FiltSong;

    [Spect, Freq, SpectTime, Power] = spectrogram(RawData, 256, 128, 256, Fs);
    FreqRows = find((Freq >= 1700) & (Freq <= 7100));

    Spect = Spect(FreqRows,:);
    Power = Power(FreqRows,:);

    SpectFs = 1/(SpectTime(2) - SpectTime(1));

    [onsets, offsets] = SegmentSong(10*log10(sum(Power)), SpectFs, 5, 10, -70);

    if (length(onsets) == 0)
        continue;
    end
    
    Notes.onsets = onsets/1000;
    Notes.offsets = offsets/1000;

    SAPFeatures = CalculateSAPFeatures(SongFileName, Fs);
    for j = 1:length(Notes.onsets),
        NoteInfo((end+1),1) = (Notes.offsets(j) - Notes.onsets(j)) * 1000;
        StartIndex = find(SAPFeatures.Time < Notes.onsets(j)/1000, 1, 'last');
        EndIndex = find(SAPFeatures.Time < Notes.offsets(j)/1000, 1, 'last');
        NoteInfo(end, 2) = mean(SAPFeatures.AM(StartIndex:EndIndex));
        NoteInfo(end, 3) = mean(SAPFeatures.Entropy(StartIndex:EndIndex));
        NoteInfo(end, 4) = mean(SAPFeatures.Freq(StartIndex:EndIndex));
        NoteInfo(end, 5) = mean(SAPFeatures.Pitch(StartIndex:EndIndex));
        NoteInfo(end, 6) = mean(SAPFeatures.FM(StartIndex:EndIndex));
        NoteInfo(end, 7) = mean(SAPFeatures.Amplitude(StartIndex:EndIndex));
        NoteInfo(end, 8) = mean(SAPFeatures.PitchGoodness(StartIndex:EndIndex));
    end
end

disp('Finished');