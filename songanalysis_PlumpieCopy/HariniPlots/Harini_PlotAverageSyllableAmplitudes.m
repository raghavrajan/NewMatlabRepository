function [] = Harini_PlotAverageSyllableAmplitudes(IndividualBirds, BirdNames, DirBout_Stats, UndirBout_Stats)

FinalFigureDir = '/home/raghav/StudentRelated/Harini/Manuscript/results';
DistanceCode = {'L0' 'L1' 'L2' 'L3' 'L4' 'UN'};

% I want to plot average syllable waveforms for L0 DIR and UN songs for a
% couple of birds specified below

BirdsToPlot = {'brwn99red65' 'brwn89red54'};
BirdTitles = {'Bird #11' 'Bird #12'};

MinTrialNo = 7;
% Now to do the log amplitude calculated similar to Kao et al.

close all;
figure;
p = panel();
p.pack({1/2 1/2});
p(1).pack('h', {1/3 1/3 1/3});
p(2).pack('h', {1/3 1/3 1/3});

for BirdIndex = 1:length(BirdsToPlot),
    i = strmatch(BirdsToPlot{BirdIndex}, BirdNames, 'exact');
    
    disp(['Bird #', num2str(i), ': ', BirdNames{i}]);
    
    % Get log amplitude waveforms for all 3 syllables for L0 dir bouts and
    % undir bouts
    for j = 1:length(DirBout_Stats{i}{1}.LogAmplitude),
        for k = 1:length(DirBout_Stats{i}{1}.LogAmplitude{j}),
            if ~isnan(DirBout_Stats{i}{1}.LogAmplitude{j}(k))
                [DataDir, SongFileName, SongFileExt] = fileparts(IndividualBirds(i).SongFileNames{IndividualBirds(i).Bouts(DirBout_Stats{i}{1}.LogAmp_BoutIndex{j}(k), 3)});
                [RawData, Fs] = GetData(DataDir, [SongFileName, SongFileExt], IndividualBirds(i).SortedBirdParameters(1).FileType, 1);
                [LogAmplitude] = ASSLCalculateLogAmplitudeKao(RawData, Fs, [], [], []);
                L0Dir_AmplitudeWaveforms{BirdIndex}{j}{k} = LogAmplitude(round(DirBout_Stats{i}{1}.LogAmp_SyllOnsets{j}(k) * Fs/1000):round(DirBout_Stats{i}{1}.LogAmp_SyllOffsets{j}(k) * Fs/1000));
            end
        end
    end

    for j = 1:length(UndirBout_Stats{i}{end}.LogAmplitude),
        for k = 1:length(UndirBout_Stats{i}{end}.LogAmplitude{j}),
            if ~isnan(UndirBout_Stats{i}{end}.LogAmplitude{j}(k))
                [DataDir, SongFileName, SongFileExt] = fileparts(IndividualBirds(i).SongFileNames{IndividualBirds(i).Bouts(UndirBout_Stats{i}{end}.LogAmp_BoutIndex{j}(k), 3)});
                [RawData, Fs] = GetData(DataDir, [SongFileName, SongFileExt], IndividualBirds(i).SortedBirdParameters(1).FileType, 1);
                [LogAmplitude] = ASSLCalculateLogAmplitudeKao(RawData, Fs, [], [], []);
                UnDir_AmplitudeWaveforms{BirdIndex}{j}{k} = LogAmplitude(round(UndirBout_Stats{i}{end}.LogAmp_SyllOnsets{j}(k) * Fs/1000):round(UndirBout_Stats{i}{end}.LogAmp_SyllOffsets{j}(k) * Fs/1000));
            end
        end
    end    
    
    % Now to warp the data to a common length
   
    for j = 1:length(L0Dir_AmplitudeWaveforms{BirdIndex}),
        L0DirLens = cellfun(@length, L0Dir_AmplitudeWaveforms{BirdIndex}{j});
        UnDirLens = cellfun(@length, UnDir_AmplitudeWaveforms{BirdIndex}{j});
        MedianLen = median([L0DirLens(:); UnDirLens(:)]);
        for k = 1:length(L0Dir_AmplitudeWaveforms{BirdIndex}{j}),
            WarpedL0Dir_AmplitudeWaveforms{BirdIndex}{j}(k,:) = interp1((1:1:length(L0Dir_AmplitudeWaveforms{BirdIndex}{j}{k}))/Fs, L0Dir_AmplitudeWaveforms{BirdIndex}{j}{k}, linspace(1/Fs, length(L0Dir_AmplitudeWaveforms{BirdIndex}{j}{k})/Fs, MedianLen));
        end
        for k = 1:length(UnDir_AmplitudeWaveforms{BirdIndex}{j}),
            WarpedUnDir_AmplitudeWaveforms{BirdIndex}{j}(k,:) = interp1((1:1:length(UnDir_AmplitudeWaveforms{BirdIndex}{j}{k}))/Fs, UnDir_AmplitudeWaveforms{BirdIndex}{j}{k}, linspace(1/Fs, length(UnDir_AmplitudeWaveforms{BirdIndex}{j}{k})/Fs, MedianLen));
        end
        p(BirdIndex,j).select();
        hold on;
        patch([(1:1:MedianLen)/Fs fliplr((1:1:MedianLen)/Fs)], [(mean(WarpedL0Dir_AmplitudeWaveforms{BirdIndex}{j}) - std(WarpedL0Dir_AmplitudeWaveforms{BirdIndex}{j})/sqrt(length(L0DirLens))) fliplr(mean(WarpedL0Dir_AmplitudeWaveforms{BirdIndex}{j}) + std(WarpedL0Dir_AmplitudeWaveforms{BirdIndex}{j})/sqrt(length(L0DirLens)))], 'k', 'FaceAlpha', 0.5, 'EdgeColor', 'none'); 
        plot((1:1:MedianLen)/Fs, mean(WarpedL0Dir_AmplitudeWaveforms{BirdIndex}{j}), 'k');
        
        patch([(1:1:MedianLen)/Fs fliplr((1:1:MedianLen)/Fs)], [(mean(WarpedUnDir_AmplitudeWaveforms{BirdIndex}{j}) - std(WarpedUnDir_AmplitudeWaveforms{BirdIndex}{j})/sqrt(length(UnDirLens))) fliplr(mean(WarpedUnDir_AmplitudeWaveforms{BirdIndex}{j}) + std(WarpedUnDir_AmplitudeWaveforms{BirdIndex}{j})/sqrt(length(UnDirLens)))], 'k', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', [0.6 0.6 0.6]); 
        plot((1:1:MedianLen)/Fs, mean(WarpedUnDir_AmplitudeWaveforms{BirdIndex}{j}), 'k', 'Color', [0.6 0.6 0.6]);
        axis tight;
        Temp = axis;
        Temp = [1/Fs MedianLen/Fs 0.94*Temp(3) 1.02*Temp(4)];
        axis(Temp);
        plot([Temp(1) Temp(1)+0.025], ones(1,2)*Temp(3), 'k', 'LineWidth', 1.5);
        set(gca, 'XColor', 'w');
        set(gca, 'YColor', 'w');
        title(['Syllable ', DirBout_Stats{i}{1}.LogAmplitude_SyllLabel(j)]);
    end        
end

p.fontsize = 12;
p.de.margin = 8;
p.marginleft = 10;
p.margintop = 10;
p.marginbottom = 5;
% p.marginright = 10;

set(gcf, 'Color', 'w');
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [10 1.4 4 2]);
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(FinalFigureDir, ['brwn99red65_brwn89red54_SyllableAverageWaveforms.eps']), '-depsc2', '-r600');
print(fullfile(FinalFigureDir, ['brwn99red65_brwn89red54_SyllableAverageWaveforms.png']), '-dpng', '-r600');

disp('Done with stats and plots');
