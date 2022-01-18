function [] = Harini_PlotRawWaveforms_LogAmplitudes(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% ======= Raw Waveforms and Log Amplitudes ===================================
% Part of the scripts used to make plots for Harini's data
% This function plots the raw waveforms and log amplitudes for individual
% syllables
% =========================================================================

Colours = 'rgbcmk';
Symbols = '+o<sd';

% First, we want to get all the log amplitudes for L0 Dir songs and
% undirected songs and find the individual syllable closest to the mean and
% also the median. Then plot the raw waveform for these syllables. Then
% plot the mean log amplitude with error bars.
for i = 1:length(IndividualBirds),
    MotifLabels = IndividualBirds(i).SortedBirdParameters(1).MotifLabels;
    for j = 1:length(MotifLabels),
        MotifLabelIndex = find(IndividualBirds(i).UniqueSyllLabels == MotifLabels(j));
        % First find all undir bouts and all L0 dir bouts
        ConditionColumnIndex = find(cellfun(@length,strfind(IndividualBirds(i).BoutStatisticsColumnNames, 'Condition')));
        UnDirectedBouts = find(IndividualBirds(i).BoutStatistics(:,ConditionColumnIndex) == 6); % 6 for undirected condition
        
        L0Bouts = find(IndividualBirds(i).BoutStatistics(:,ConditionColumnIndex) == 1); % 1 for L0 condition
        DirectedL0Bouts = intersect(find(cellfun(@length, strfind(IndividualBirds(i).BoutCategorisation, 'D'))), L0Bouts);
        
        L0DirAmplitudes = [];
        L0DirBoutIndices = [];
        for k = DirectedL0Bouts(:)',
            if (~isempty(IndividualBirds(i).AllSylls(MotifLabelIndex).LogAmplitude{k}))
                L0DirAmplitudes = [L0DirAmplitudes; (IndividualBirds(i).AllSylls(MotifLabelIndex).LogAmplitude{k}(:))];
                L0DirBoutIndices = [L0DirBoutIndices; [ones(size((IndividualBirds(i).AllSylls(MotifLabelIndex).LogAmplitude{k}(:))))*k (1:1:length((IndividualBirds(i).AllSylls(MotifLabelIndex).LogAmplitude{k}(:))))']];
            end
        end
        % Eliminate nans here
        NanIndices = find(isnan(L0DirAmplitudes));
        L0DirAmplitudes(NanIndices) = [];
        L0DirBoutIndices(NanIndices,:) = [];

        
        UnDirectedAmplitudes = [];
        UnDirectedBoutIndices = [];
        for k = UnDirectedBouts(:)',
            if (~isempty(IndividualBirds(i).AllSylls(MotifLabelIndex).LogAmplitude{k}))
                UnDirectedAmplitudes = [UnDirectedAmplitudes; (IndividualBirds(i).AllSylls(MotifLabelIndex).LogAmplitude{k}(:))];
                UnDirectedBoutIndices = [UnDirectedBoutIndices; [ones(size((IndividualBirds(i).AllSylls(MotifLabelIndex).LogAmplitude{k}(:))))*k (1:1:length((IndividualBirds(i).AllSylls(MotifLabelIndex).LogAmplitude{k}(:))))']];
            end
        end
        
        % Now that we have got all the amplitudes, have to first check
        % whether there are both L0 and undir renditions for this syllable
        if ((length(L0DirAmplitudes) < 5) || (length(UnDirectedAmplitudes) < 5))
            disp([BirdNames{i}, ': not enough data; # of dir bouts = ', num2str(length(L0DirAmplitudes)), ' and # of undir bouts = ', num2str(length(UnDirectedAmplitudes))]);
            continue;
        end
        
        [h, p] = ttest2(L0DirAmplitudes(find(~isnan(L0DirAmplitudes))), UnDirectedAmplitudes(find(~isnan(UnDirectedAmplitudes))));
        
        % Now we need to get file nos. and onsets and offsets for the
        % individual syllables
        BoutIndexColumnIndex = find(cellfun(@length,strfind(IndividualBirds(i).BoutStatisticsColumnNames, 'BoutIndex')));
        clear L0DirSyllFileName L0DirSyllOnset L0DirSyllOffset L0DirSyllLogAmplitude L0DirActualRawSyllable;
        for k = 1:length(L0DirAmplitudes),
            BoutIndex = IndividualBirds(i).BoutStatistics(L0DirBoutIndices(k,1), BoutIndexColumnIndex);
            FileIndex = IndividualBirds(i).Bouts(BoutIndex,3);
            L0DirSyllFileName{k} = IndividualBirds(i).SongFileNames{FileIndex};
            
            BoutStartSyllIndex = IndividualBirds(i).Bouts(BoutIndex,1);
            BoutEndSyllIndex = IndividualBirds(i).Bouts(BoutIndex,2);
            
            BoutSyllLabels = char(IndividualBirds(i).AllSyllableData(BoutStartSyllIndex:BoutEndSyllIndex,1));
            BoutSyllOnsets = IndividualBirds(i).AllSyllableData(BoutStartSyllIndex:BoutEndSyllIndex,4);
            BoutSyllOffsets = IndividualBirds(i).AllSyllableData(BoutStartSyllIndex:BoutEndSyllIndex,5);
            
            Matches = find(BoutSyllLabels == MotifLabels(j));
            L0DirSyllOnsetsOffsets(k,:) = [BoutSyllOnsets(Matches(L0DirBoutIndices(k,2))) BoutSyllOffsets(Matches(L0DirBoutIndices(k,2)))];
            
            [DataDir, Name, Ext] = fileparts(L0DirSyllFileName{k});
            [RawData, Fs] = GetData(DataDir, [Name, Ext], IndividualBirds(i).SortedBirdParameters(1).FileType, 0);
            
            Padding = 200; % Take 200ms of extra data
            L0DirActualRawSyllable{k} = RawData(round((L0DirSyllOnsetsOffsets(k,1)) * Fs/1000):round((L0DirSyllOnsetsOffsets(k,2)) * Fs/1000));
            RawData = RawData(round((L0DirSyllOnsetsOffsets(k,1) - Padding) * Fs/1000):round((L0DirSyllOnsetsOffsets(k,2) + Padding) * Fs/1000));
            
            [Feats, RawFeats, FeatsFs] = ASSLCalculateSAPFeatsWithOnsets(RawData, (1:1:length(RawData))/Fs, Fs, Padding/1000, Padding/1000 + (L0DirSyllOnsetsOffsets(k,2) - L0DirSyllOnsetsOffsets(k,1))/1000);
            L0DirSyllLogAmplitude{k} = RawFeats.LogAmplitude{1};
        end
        
        [SortedVals, SortedIndices] = sort(L0DirAmplitudes);
        L0DirAmplitudes = L0DirAmplitudes(SortedIndices);
        L0DirActualRawSyllable = L0DirActualRawSyllable(SortedIndices);
        L0DirSyllLogAmplitude = L0DirSyllLogAmplitude(SortedIndices);
        
        % Now we need to get file nos. and onsets and offsets for the
        % individual syllables
        BoutIndexColumnIndex = find(cellfun(@length,strfind(IndividualBirds(i).BoutStatisticsColumnNames, 'BoutIndex')));
        clear UnDirectedSyllFileName UnDirectedSyllOnset UnDirectedSyllOffset UnDirectedSyllLogAmplitude UnDirectedActualRawSyllable;
        for k = 1:length(UnDirectedAmplitudes),
            BoutIndex = IndividualBirds(i).BoutStatistics(UnDirectedBoutIndices(k,1), BoutIndexColumnIndex);
            FileIndex = IndividualBirds(i).Bouts(BoutIndex,3);
            UnDirectedSyllFileName{k} = IndividualBirds(i).SongFileNames{FileIndex};
            
            BoutStartSyllIndex = IndividualBirds(i).Bouts(BoutIndex,1);
            BoutEndSyllIndex = IndividualBirds(i).Bouts(BoutIndex,2);
            
            BoutSyllLabels = char(IndividualBirds(i).AllSyllableData(BoutStartSyllIndex:BoutEndSyllIndex,1));
            BoutSyllOnsets = IndividualBirds(i).AllSyllableData(BoutStartSyllIndex:BoutEndSyllIndex,4);
            BoutSyllOffsets = IndividualBirds(i).AllSyllableData(BoutStartSyllIndex:BoutEndSyllIndex,5);
            
            Matches = find(BoutSyllLabels == MotifLabels(j));
            UnDirectedSyllOnsetsOffsets(k,:) = [BoutSyllOnsets(Matches(UnDirectedBoutIndices(k,2))) BoutSyllOffsets(Matches(UnDirectedBoutIndices(k,2)))];
            
            [DataDir, Name, Ext] = fileparts(UnDirectedSyllFileName{k});
            [RawData, Fs] = GetData(DataDir, [Name, Ext], IndividualBirds(i).SortedBirdParameters(1).FileType, 0);
            
            Padding = 200; % Take 200ms of extra data
            UnDirectedActualRawSyllable{k} = RawData(round((UnDirectedSyllOnsetsOffsets(k,1)) * Fs/1000):round((UnDirectedSyllOnsetsOffsets(k,2)) * Fs/1000));
            RawData = RawData(round((UnDirectedSyllOnsetsOffsets(k,1) - Padding) * Fs/1000):round((UnDirectedSyllOnsetsOffsets(k,2) + Padding) * Fs/1000));
            
            [Feats, RawFeats, FeatsFs] = ASSLCalculateSAPFeatsWithOnsets(RawData, (1:1:length(RawData))/Fs, Fs, Padding/1000, Padding/1000 + (UnDirectedSyllOnsetsOffsets(k,2) - UnDirectedSyllOnsetsOffsets(k,1))/1000);
            UnDirectedSyllLogAmplitude{k} = RawFeats.LogAmplitude{1};
        end
        [SortedVals, SortedIndices] = sort(UnDirectedAmplitudes);
        UnDirectedAmplitudes = UnDirectedAmplitudes(SortedIndices);
        UnDirectedActualRawSyllable = UnDirectedActualRawSyllable(SortedIndices);
        UnDirectedSyllLogAmplitude = UnDirectedSyllLogAmplitude(SortedIndices);
       
        figure;
        set(gcf, 'Color', 'w');
        set(gcf, 'Position', [680 90 600 900]);
        subplot(2,1,1);
        hold on;
        % First plot raw waveforms
        L0DirMedianIndex = round(length(find(~isnan(L0DirAmplitudes)))/2);
        UnDirMedianIndex = round(length(find(~isnan(UnDirectedAmplitudes)))/2);
        plot((1:1:length(L0DirActualRawSyllable{L0DirMedianIndex}))/Fs, L0DirActualRawSyllable{L0DirMedianIndex}, 'r');
        plot((length(L0DirActualRawSyllable{L0DirMedianIndex})/Fs)*1.05 + (1:1:length(UnDirectedActualRawSyllable{UnDirMedianIndex}))/Fs, UnDirectedActualRawSyllable{UnDirMedianIndex}, 'b');
        axis tight;
        set(gca, 'FontSize', 14);
        xlabel('Time (msec)');
        ylabel('Amplitude (V)');
        if (nanmean(L0DirAmplitudes) > nanmean(UnDirectedAmplitudes))
            if (p < 0.05)
                title([BirdNames{i}, ': Syll ', MotifLabels(j), ': Dir > Undir *']);
            else
                title([BirdNames{i}, ': Syll ', MotifLabels(j), ': Undir > Dir']);
            end
        else
            if (p < 0.05)
                title([BirdNames{i}, ': Syll ', MotifLabels(j), ': Undir > Dir *']);
            else
                title([BirdNames{i}, ': Syll ', MotifLabels(j), ': Dir > Undir']);
            end
        end
        
        subplot(2,1,2);
        hold on;
        % Now warp the data and plot the average and sem of the waveforms
        MedianLen = median([cellfun(@length, L0DirSyllLogAmplitude(find(~isnan(L0DirAmplitudes)))) cellfun(@length, UnDirectedSyllLogAmplitude(find(~isnan(UnDirectedAmplitudes))))]);
        MedianTime = median([cellfun(@length, L0DirActualRawSyllable(find(~isnan(L0DirAmplitudes)))) cellfun(@length, UnDirectedActualRawSyllable(find(~isnan(UnDirectedAmplitudes))))]);
        L0DirWarpedAmplitudes = [];
        UnDirectedWarpedAmplitudes = [];
        for k = 1:length(L0DirAmplitudes),
            if (~isnan(L0DirAmplitudes(k)))
                L0DirWarpedAmplitudes(k,:) = interp1((1:1:length(L0DirSyllLogAmplitude{k}))/Fs, L0DirSyllLogAmplitude{k}, linspace(1/Fs, length(L0DirSyllLogAmplitude{k})/Fs, MedianLen));
            end
        end
        for k = 1:length(UnDirectedAmplitudes),
            if (~isnan(UnDirectedAmplitudes(k)))
                UnDirectedWarpedAmplitudes(k,:) = interp1((1:1:length(UnDirectedSyllLogAmplitude{k}))/Fs, UnDirectedSyllLogAmplitude{k}, linspace(1/Fs, length(UnDirectedSyllLogAmplitude{k})/Fs, MedianLen));
            end
        end
        Time = 1000 * linspace(1/Fs, MedianTime/Fs, MedianLen);
        patch([Time fliplr(Time)], [(mean(L0DirWarpedAmplitudes) + std(L0DirWarpedAmplitudes)/sqrt(size(L0DirWarpedAmplitudes, 1))) fliplr(mean(L0DirWarpedAmplitudes) - std(L0DirWarpedAmplitudes)/sqrt(size(L0DirWarpedAmplitudes, 1)))], 'k', 'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.5);
        plot(Time, mean(L0DirWarpedAmplitudes), 'r');
        
        patch([Time fliplr(Time)], [(mean(UnDirectedWarpedAmplitudes) + std(UnDirectedWarpedAmplitudes)/sqrt(size(UnDirectedWarpedAmplitudes, 1))) fliplr(mean(UnDirectedWarpedAmplitudes) - std(UnDirectedWarpedAmplitudes)/sqrt(size(UnDirectedWarpedAmplitudes, 1)))], 'k', 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.5);
        plot(Time, mean(UnDirectedWarpedAmplitudes), 'b');
        axis tight;
        set(gca, 'FontSize', 14);
        xlabel('Time (msec)');
        ylabel('Log amplitude (dB)');
    end
end

disp('Finished');