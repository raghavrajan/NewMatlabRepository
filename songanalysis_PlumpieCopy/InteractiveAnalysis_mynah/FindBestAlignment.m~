function [] = FindBestAlignment(DirFileInfo, UnDirFileInfo, MedianMotif, MotifString, DirectoryName, FileType, SongChannelNo)

[RawData, Fs] = ReadOKrankData(DirectoryName, MedianMotif.FileName{1}, SongChannelNo);
[b, a] = butter(8, [300*2/Fs, 10000*2/Fs]);
RawData = filtfilt(b, a, RawData);

for i = 1:length(MedianMotif.SyllableLengths),
    if (i == 1)
        StartIndex = round((MedianMotif.StartTime + MedianMotif.SyllableStartings(i) - 0.04) * Fs);
    else
        StartIndex = round((MedianMotif.StartTime + MedianMotif.SyllableStartings(i) - (MedianMotif.GapLengths(i-1)/2)) * Fs);
    end
    if (i == length(MedianMotif.SyllableLengths))
        EndIndex = round((MedianMotif.StartTime + MedianMotif.SyllableStartings(i) + MedianMotif.SyllableLengths(i) + 0.04) * Fs);
    else
        EndIndex = round((MedianMotif.StartTime + MedianMotif.SyllableStartings(i) + MedianMotif.SyllableLengths(i) + (MedianMotif.GapLengths(i)/2)) * Fs);
    end
    SyllableData = RawData(StartIndex:EndIndex);
    [Spect, Freq, Time] = spectrogram(SyllableData, 128, 124, 128, Fs);
    FreqRows = find((Freq >= 1700) & (Freq <= 7100));
    Spect = Spect(FreqRows,:);
    MeanSpect = (mean(mean(Spect)));
    STDSpect = (sqrt((sum(sum((Spect - MeanSpect).*(Spect - MeanSpect))))/((size(Spect, 1) * size(Spect, 2)) - 1)));
    Spect = (Spect - mean(mean(Spect)))/STDSpect;
    
    for j = 2:size(Spect, 2),
        SpectDiff(:,(j-1)) = Spect(:,j) - Spect(:,(j-1));
    end
    for j = 1:size(DirFileInfo.Syllables.Length, 1),
        [RawData2, Fs2] = ReadOKrankData(DirectoryName, DirFileInfo.FileNames{DirFileInfo.Syllables.Index(j)}, SongChannelNo);
        [b, a] = butter(8, [300*2/Fs2, 10000*2/Fs2]);
        RawData2 = filtfilt(b, a, RawData2);
        if (i == 1)
            StartIndex = round((DirFileInfo.Syllables.Start(j,i) - 0.04) * Fs);
        else
            StartIndex = round((DirFileInfo.Syllables.Start(j,i) - (DirFileInfo.Gaps.Length(j,(i-1))/2)) * Fs);
        end
        if (i == length(MedianMotif.SyllableLengths))
            EndIndex = round((DirFileInfo.Syllables.End(j,i) + 0.04) * Fs);
        else
            EndIndex = round((DirFileInfo.Syllables.End(j,i) + (DirFileInfo.Gaps.Length(j,i)/2)) * Fs);
        end
        
        SyllableData2 = RawData2(StartIndex:EndIndex);
        [Spect2, Freq2, Time2] = spectrogram(SyllableData2, 128, 124, 128, Fs2);
        FreqRows2 = find((Freq2 >= 1700) & (Freq2 <= 7100));
        Spect2 = Spect2(FreqRows2,:);
        MeanSpect2 = (mean(mean(Spect2)));
        STDSpect2 = (sqrt((sum(sum((Spect2 - MeanSpect2).*(Spect2 - MeanSpect2))))/((size(Spect2, 1) * size(Spect2, 2)) - 1)));
        Spect2 = (Spect2 - mean(mean(Spect2)))/STDSpect2;
        
        clear SpectDiff2;
        for k = 2:size(Spect2, 2),
            SpectDiff2(:,(k-1)) = Spect2(:,k) - Spect2(:,(k-1));
        end
    end
end

disp('Finding best alignment');
