function [AronovFeeLogAmpValues, BaselineAmpValue, SyllOnsetValues, SyllOffsetValues, FirstHalfAmpValue, SecondHalfAmpValue] = LSINA_CalcAronovFeeLogAmpValues(BirdParameters)

if (BirdParameters.Continuousdata == 0)
    for i = 1:length(BirdParameters.SongFileNames),
        if (isempty(BirdParameters.NoteInfo{i}.onsets))
            AronovFeeLogAmpValues{i} = [];
            BaselineAmpValue{i} = [];
            SyllOnsetValues{i} = [];
            SyllOffsetValues{i} = [];
            FirstHalfAmpValue{i} = [];
            SecondHalfAmpValue{i} = [];
            continue;
        end
        [RawData, Fs] = ASSLGetRawData(BirdParameters.DataDirectory, BirdParameters.SongFileNames{i}, BirdParameters.FileType, 0);

        Time = (1:1:length(RawData))/Fs;
        FFTWinSize = 5; % in ms
        F_High = 1000;
        F_Low = 8000;
        [LogAmplitude] = ASSLCalculateLogAmplitudeAronovFee_WithoutLog(RawData, Fs, Time, FFTWinSize, [], F_High, F_Low);
        Temp = gmdistribution.fit(LogAmplitude(:), 2);
        for j = 1:length(BirdParameters.NoteInfo{i}.onsets),
            SyllableOnsetTimeIndex = round(BirdParameters.NoteInfo{i}.onsets(j)/1000 * Fs);
            SyllableOffsetTimeIndex = round(BirdParameters.NoteInfo{i}.offsets(j)/1000 * Fs);
            if (SyllableOnsetTimeIndex(1) <= 0)
                SyllableOnsetTimeIndex(1) = 1;
            end
            if (SyllableOffsetTimeIndex > length(RawData))
                SyllableOffsetTimeIndex = length(RawData);
            end
            AronovFeeLogAmpValues{i}(j) = mean(LogAmplitude(SyllableOnsetTimeIndex:SyllableOffsetTimeIndex));
            BaselineAmpValue{i}(j) = min(Temp.mu);
            HalfwayIndex = round(mean([SyllableOnsetTimeIndex SyllableOffsetTimeIndex]));
            FirstHalfAmpValue{i}(j) = mean(LogAmplitude(SyllableOnsetTimeIndex:HalfwayIndex));
            SecondHalfAmpValue{i}(j) = mean(LogAmplitude(HalfwayIndex:SyllableOffsetTimeIndex));
            SyllOnsetValues{i}(j) = LogAmplitude(SyllableOnsetTimeIndex);
            SyllOffsetValues{i}(j) = LogAmplitude(SyllableOffsetTimeIndex);
        end
    end
end
