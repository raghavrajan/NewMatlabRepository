function [] = PrepTemplates(NoteFileList, RawDataDir, SongChanNo, FileType, OutputFileName)

Labels = [];
for i = 1:length(NoteFileList),
    Notes{i} = load(NoteFileList(i).name);
    Labels = unique([Labels Notes{i}.labels]);
    if (strfind(FileType, 'okrank'))
        [Song{i}, Fs{i}] = ReadOKrankData(RawDataDir, NoteFileList(i).name(1:(end-8)), SongChanNo);
    end
end

FreqVect = [300:75:10000];
for i = 1:length(Labels),
    TempPSD = zeros(length(FreqVect), 8);
    TempDurations = [];
    Index = 0;
    for j = 1:length(NoteFileList),
        LabelIndices = find(Notes{j}.labels == Labels(i));
        for k = 1:length(LabelIndices),
            Time = ((1:1:length(Song{j}))/Fs{j}) * 1000;
            OnsetIndices = find(Time <= Notes{j}.onsets(LabelIndices(k)), 1, 'last');
            OffsetIndices = find(Time <= Notes{j}.offsets(LabelIndices(k)), 1, 'last');
            [S, F, T, P] = spectrogram(Song{j}(OnsetIndices:OffsetIndices),[], 0, FreqVect, Fs{j});
            TempPSD = TempPSD + P;
            Index = Index + 1;
            TempDurations = [TempDurations; ((Notes{j}.offsets(LabelIndices(k)) - Notes{j}.onsets(LabelIndices(k))) / 1000)];
       end
    end
    if (Index > 0)
        TempPSD = TempPSD/Index;
        MeanTempPSD = mean(mean(TempPSD));
        STDTempPSD = sqrt(sum(sum((TempPSD - MeanTempPSD).*(TempPSD - MeanTempPSD)))/((size(TempPSD,1) * size(TempPSD, 2)) - 1));
        TempPSD = (TempPSD - MeanTempPSD)/STDTempPSD;
        PSD{i} = TempPSD;
        Durations{i} = mean(TempDurations);
    end
end

save(OutputFileName, 'PSD', 'Labels', 'NoteFileList', 'Durations');
disp('Finished');