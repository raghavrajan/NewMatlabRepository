function [] = MakeSyllDistributionHistogramsAronovFee(RawDataDirs, FileLists, FileType, SongChanNo, OutputFileName, BirdName, Legend)

PresentDir = pwd;

if (exist(OutputFileName, 'file'))
    load(OutputFileName);
else
    for i = 1:length(FileLists),
        cd(PresentDir);
        Fid = fopen(FileLists{i}, 'r');
        Temp = textscan(Fid, '%s', 'DeLimiter', '\n');
        SongFiles{i} = Temp{1};
        fclose(Fid);
        
        Durs{i} = [];
        ConsecutiveSyllDurs{i} = [];
        ConsecutiveGaps{i} = [];
        
        NextGap{i} = [];
        PrevGap{i} = [];
        
        for j = 1:length(SongFiles{i}),
            [RawData, Fs] = ASSLGetRawData(RawDataDirs{i}, SongFiles{i}{j}, FileType, SongChanNo);

            Time = (1:1:length(RawData))/Fs;

            [LogAmplitude] = ASSLCalculateLogAmplitudeAronovFee(RawData, Fs, Time, 8, 0.5);

            Threshold{i}{j} = ASSLCalculateFisherThreshold(LogAmplitude);
            [SyllOnsets{i}{j}, SyllOffsets{i}{j}] = ASSLSegmentDataAronovFee(LogAmplitude, Fs, 7, 7, Threshold{i}{j});
            TempDurs = SyllOffsets{i}{j} - SyllOnsets{i}{j};
            TempGaps = SyllOnsets{i}{j}(2:end) - SyllOffsets{i}{j}(1:end-1);
            TempDurs = TempDurs(:);
            TempGaps = TempGaps(:);
            
            Durs{i} = [Durs{i}; TempDurs(:)];
            ConsecutiveSyllDurs{i} = [ConsecutiveSyllDurs{i}; [TempDurs(1:end-1) TempDurs(2:end)]];
            NextGap{i} = [NextGap{i}; [TempDurs(1:end-1) TempGaps]];
            PrevGap{i} = [PrevGap{i}; [TempDurs(2:end) TempGaps]];
            
            ConsecutiveGaps{i} = [ConsecutiveGaps{i}; [TempGaps(1:end-1) TempGaps(2:end)]];
            disp([' Detected ', num2str(length(SyllOnsets{i}{j})), ' syllables for ', SongFiles{i}{j}]);
        end
    end
end

cd(PresentDir);
save(OutputFileName);

% First plot syllable duration histograms

figure;
set(gcf, 'Position', [418 222 800 450]);
hold on;
Colours = 'rgbcmy';
for i = 1:length(Durs),
    plot([0:5:max(Durs{i})], histc(Durs{i}, [0:5:max(Durs{i})]) * 100/sum(histc(Durs{i}, [0:5:max(Durs{i})])), Colours(i));
end

set(gcf, 'Color', 'k');
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
xlabel('Duration (ms)', 'FontSize', 14, 'FontWeight', 'normal', 'FontName', 'Arial');
ylabel('%', 'FontSize', 14, 'FontWeight', 'normal', 'FontName', 'Arial');
title([BirdName, ': Undir syllable duration histograms'], 'FontSize', 16, 'Color', 'w');
LegendHandle = legend(Legend);
set(LegendHandle, 'TextColor', 'w', 'FontSize', 14);

% Next plot the syllable duration contingencies
figure;
set(gcf, 'Position', [418 222 800 450]);
hold on;
Colours = 'rgbcmy';

for i = 1:length(ConsecutiveSyllDurs),
    Edges = 0:25:500;
    for j = 1:length(Edges)-1,
        Indices = find((ConsecutiveSyllDurs{i}(:,1) >= Edges(j)) & (ConsecutiveSyllDurs{i}(:,1) < Edges(j+1)));
        SyllDurContingencies{i}(:,j) = histc(ConsecutiveSyllDurs{i}(Indices,2), Edges)/sum(histc(ConsecutiveSyllDurs{i}(Indices,2), Edges));
        
        Indices = find((NextGap{i}(:,1) >= Edges(j)) & (NextGap{i}(:,1) < Edges(j+1)));
        NextGapContingencies{i}(:,j) = histc(NextGap{i}(Indices,2), Edges)/sum(histc(NextGap{i}(Indices,2), Edges));
        
        Indices = find((ConsecutiveGaps{i}(:,1) >= Edges(j)) & (ConsecutiveGaps{i}(:,1) < Edges(j+1)));
        GapContingencies{i}(:,j) = histc(ConsecutiveGaps{i}(Indices,2), Edges)/sum(histc(ConsecutiveGaps{i}(Indices,2), Edges));
        
        Indices = find((PrevGap{i}(:,1) >= Edges(j)) & (PrevGap{i}(:,1) < Edges(j+1)));
        PrevGapContingencies{i}(:,j) = histc(PrevGap{i}(Indices,2), Edges)/sum(histc(PrevGap{i}(Indices,2), Edges));
    end
end
    
set(gcf, 'Color', 'k');
set(gca, 'Color', 'k');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
xlabel('Duration (ms)', 'FontSize', 14, 'FontWeight', 'normal', 'FontName', 'Arial');
ylabel('%', 'FontSize', 14, 'FontWeight', 'normal', 'FontName', 'Arial');
title([BirdName, ': Undir syllable duration histograms'], 'FontSize', 16, 'Color', 'w');
LegendHandle = legend(Legend);
set(LegendHandle, 'TextColor', 'w', 'FontSize', 14);