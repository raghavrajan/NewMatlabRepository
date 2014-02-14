function [SyllData] = ASSLGetIndividualSyllData(handles)

PreTime = 25;
PostTime = 25;

for i = 1:length(handles.ASB.UniqueSyllLabels),
    Matches = find(handles.ASB.AllSyllLabels == handles.ASB.UniqueSyllLabels{i});
    SyllData{i}.Label = handles.ASB.UniqueSyllLabels{i};
    for j = 1:length(Matches),
        SyllData{i}.FileSyllIndices(j,:) = [handles.ASB.AllSyllIndices(Matches(j),1) handles.ASB.AllSyllIndices(Matches(j),2)];
        if (j > 1)
            if (handles.ASB.AllSyllIndices(Matches(j),1) ~= handles.ASB.AllSyllIndices(Matches(j-1),1))
                [RawData, Fs] = ASSLGetRawData(handles.ASB.DirName, handles.ASB.FileName{handles.ASB.AllSyllIndices(Matches(j),1)}, handles.ASB.FileType, handles.ASB.SongChanNo);
                Time = (1:1:length(RawData))/Fs * 1000;
                [LogAmplitude] = ASSLCalculateLogAmplitude(RawData, Fs, Time, handles.ASB.FFTWinSizeSegmenting, handles.ASB.FFTWinOverlapSegmenting);
            end
        else
            [RawData, Fs] = ASSLGetRawData(handles.ASB.DirName, handles.ASB.FileName{handles.ASB.AllSyllIndices(Matches(j),1)}, handles.ASB.FileType, handles.ASB.SongChanNo);
            Time = (1:1:length(RawData))/Fs * 1000;
            [LogAmplitude] = ASSLCalculateLogAmplitude(RawData, Fs, Time, handles.ASB.FFTWinSizeSegmenting, handles.ASB.FFTWinOverlapSegmenting);
        end
        
        StartIndex = find(Time >= (handles.ASB.SyllOnsets{handles.ASB.AllSyllIndices(Matches(j),1)}(handles.ASB.AllSyllIndices(Matches(j),2)) - PreTime), 1, 'first');
        EndIndex = find(Time <= (handles.ASB.SyllOffsets{handles.ASB.AllSyllIndices(Matches(j),1)}(handles.ASB.AllSyllIndices(Matches(j),2)) + PostTime), 1, 'last');
        SyllData{i}.Amplitude{j} = LogAmplitude(StartIndex:EndIndex);
        SyllData{i}.Time{j} = Time(StartIndex:EndIndex) - handles.ASB.SyllOnsets{handles.ASB.AllSyllIndices(Matches(j),1)}(handles.ASB.AllSyllIndices(Matches(j),2));
        SyllData{i}.OffsetTime{j} = SyllData{i}.Time{j}(end) - PostTime;
        SyllData{i}.OffsetAlignedTime{j} = Time(StartIndex:EndIndex) - handles.ASB.SyllOffsets{handles.ASB.AllSyllIndices(Matches(j),1)}(handles.ASB.AllSyllIndices(Matches(j),2));
    end
end

disp('Finished getting individual syllable data');