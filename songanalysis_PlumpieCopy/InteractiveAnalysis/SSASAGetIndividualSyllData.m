function [DataStruct] = SSASAGetIndividualSyllData(InputData)

PreTime = InputData.SSASA.PreTime;
PostTime = InputData.SSASA.PostTime;

for i = 1:length(InputData.SSASA.SyllLabels),
    DataStruct{i}.SyllLabel = InputData.SSASA.SyllLabels{i};
    SyllIndex = 0;
    for j = 1:length(InputData.SSASA.UnDirFileInfo.NoteLabels),
        Matches = find(InputData.SSASA.UnDirFileInfo.NoteLabels{j} == DataStruct{i}.SyllLabel);
        if (~isempty(Matches))
            for k = 1:length(Matches),
                SyllIndex = SyllIndex + 1;
                DataStruct{i}.SyllStartTime{SyllIndex} = InputData.SSASA.UnDirFileInfo.NoteOnsets{j}(Matches(k));
                DataStruct{i}.SyllEndTime{SyllIndex} = InputData.SSASA.UnDirFileInfo.NoteOffsets{j}(Matches(k));
                DataStruct{i}.SyllDur{SyllIndex} = [InputData.SSASA.UnDirFileInfo.NoteOnsets{j}(Matches(k)) InputData.SSASA.UnDirFileInfo.NoteOffsets{j}(Matches(k))];
                DataStruct{i}.SyllDuration{SyllIndex} = DataStruct{i}.SyllEndTime{SyllIndex} - DataStruct{i}.SyllStartTime{SyllIndex};
                DataStruct{i}.FileName{SyllIndex} = InputData.SSASA.UnDirFileInfo.FileNames{j};
                DataStruct{i}.FileDur{SyllIndex} = [0 InputData.SSASA.UnDirFileInfo.RecordLengths{j}];
                DataStruct{i}.SpikeTimes{SyllIndex} = InputData.SSASA.UnDirFileInfo.SpikeData.Times{j};
            end
        end
    end
    MedianSyllDuration = median(cell2mat(DataStruct{i}.SyllDuration));
    [SortedVals, SortedIndices] = sort(cell2mat(DataStruct{i}.SyllDuration));
    UnWarpedRaster = [];
    WarpedRaster = [];
    Index = 0;
    for j = SortedIndices,
        TempSpikeTimes = DataStruct{i}.SpikeTimes{j}(find((DataStruct{i}.SpikeTimes{j} >= (DataStruct{i}.SyllStartTime{j} - PreTime)) & (DataStruct{i}.SpikeTimes{j} < (DataStruct{i}.SyllEndTime{j} + PostTime))));
        TempSpikeTimes = TempSpikeTimes(:);
        TempSpikeTimes = TempSpikeTimes - DataStruct{i}.SyllStartTime{j};
        UnWarpedRaster = [UnWarpedRaster; [TempSpikeTimes ones(size(TempSpikeTimes))*Index]];
        Index = Index + 1;
    end
    DataStruct{i}.UnWarpedRaster = UnWarpedRaster;
    DataStruct{i}.UnWarpedRasterAxis = [-PreTime (DataStruct{i}.SyllDuration{j} + PostTime) -1 (Index + 1)];
end
disp('Finished parsing data into separate syllables');