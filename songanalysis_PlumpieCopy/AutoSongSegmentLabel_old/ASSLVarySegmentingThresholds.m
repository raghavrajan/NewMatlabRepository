function [] = ASSLVarySegmentingThresholds(FileList, Threshold, MinInt, MinDur, SmoothWindow, WindowOverlap, FileType, SongChanNo)

Fid = fopen(FileList, 'r');
Temp = textscan(Fid, '%s', 'DeLimiter', '\n');
Temp = Temp{1}(1:50);
fclose(Fid);

for i = 1:length(Temp),
    SlashIndex = find((Temp{i} == '/') | (Temp{i} == '\'));
    if (~isempty(SlashIndex))
        Temp{i} = Temp{i}(SlashIndex(end)+1:end);
    end
end


for i = 1:length(Threshold),
    SyllDur{i} = [];
    IntervalDur{i} = [];
    SAPFeats{i} = [];
end

for j = 1:length(Temp),
    [RawData, Fs] = ASSLGetRawData([pwd,'/'], Temp{j}, FileType, SongChanNo);
    Time = (1:1:length(RawData))*1000/Fs;
    LogAmplitude = ASSLCalculateLogAmplitude(RawData, Fs, [1:1:length(RawData)]/Fs, SmoothWindow, WindowOverlap);
    for i = 1:length(Threshold),
        [Onsets, Offsets] = ASSLSegmentData(LogAmplitude, Fs, MinInt, MinDur, Threshold(i));
        NumSyll{j}(i) = length(Onsets);
        if (~isempty(Onsets))
            SyllDur{i} = [SyllDur{i}; [Offsets - Onsets]];
            IntervalDur{i} = [IntervalDur{i}; [Onsets(2:end) - Offsets(1:end-1)]];
            %TempFeats = ASSLCalculateSAPFeatsWithOnsets(RawData, [1:1:length(RawData)]*1000/Fs, Fs, Onsets, Offsets);
            TempFeats = [];
            for k = 1:length(Onsets),
                StartIndex = find(Time >= Onsets(k), 1, 'first');
                EndIndex = find(Time <= Offsets(k), 1, 'last');
                TempFeats.Duration(k) = Offsets(k) - Onsets(k);
                TempFeats.LogAmplitude(k) = mean(LogAmplitude(StartIndex:EndIndex));
            end
            SAPFeats{i} = [SAPFeats{i}; TempFeats];
        end
    end
end

figure;
for i = 1:length(Threshold),
    subplot(4,4,i);
    if (~isempty([SAPFeats{i}.Duration]))
        plot([SAPFeats{i}.Duration], [SAPFeats{i}.LogAmplitude], 'k+');
        title(Threshold(i));
        axis([0 500 20 80]);
    end
end
disp('Finished');