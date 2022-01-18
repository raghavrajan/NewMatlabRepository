function [] = FindBoutStartBoutEndIntervals(InterBoutInterval, SavedDataDir)


% File name of the file used to get the data from
SongDetailsFile = '/home/raghav/StudentRelated/Harini/Manuscript/JEBManuscript/ScriptsForFigures/ContinuousDataSongAnalysis_BirdDetails_ForHariniDataRelatedAnalysis.csv';
OutputDir = '/home/raghav/StudentRelated/Harini/Manuscript/JEBManuscript/Figures';

%% Get all data
[BirdParameters, Flag] = ProcessSongData_IntoBouts(SongDetailsFile, InterBoutInterval, SavedDataDir);


%% Now, to find the fraction of intervals at the beginning greater than different lengths
Edges = 0:50:InterBoutInterval; % In ms
for i = 1:length(BirdParameters),
    ValidSongBouts = find((BirdParameters(i).Bouts(:,7) == 1) & (BirdParameters(i).Bouts(:,8) > 0) & (BirdParameters(i).Bouts(:,9) > 1));
    for j = 1:length(ValidSongBouts),
        FirstMotifSyll = NaN;
        for k = 1:length(BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Labels),
            if (~isempty(find(BirdParameters(i).MotifLabels == BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Labels(k))))
                FirstMotifSyll = k;
                break;
            end
        end
        if (isempty(find(isnan(FirstMotifSyll))))
            IntervalsBeforeFirstMotifSyll{i}{j} = BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Onsets(2:FirstMotifSyll) - BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Offsets(1:FirstMotifSyll-1);
            BinnedIntervalsBeforeFirstMotifSyll{i}(j,:) = histc(BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Onsets(2:FirstMotifSyll) - BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Offsets(1:FirstMotifSyll-1), Edges);
        end
    end
end

%% Now to find the fraction of bouts that have atleast 1 interval > each of the values in edges - sort of a cumulative fraction
for i = 1:length(BinnedIntervalsBeforeFirstMotifSyll),
    TempIntervalsBeforeFirstMotifSyll = BinnedIntervalsBeforeFirstMotifSyll{i};
    % Threshold it to 1 and 0 for anything that is > 1
    TempIntervalsBeforeFirstMotifSyll = TempIntervalsBeforeFirstMotifSyll > 0;
    % Flip lr to get fraction of intervals > a particular value
    TempIntervalsBeforeFirstMotifSyll = fliplr(TempIntervalsBeforeFirstMotifSyll);
    % Now do a cumlative sum along the row and threshold it and divide by
    % total number of bouts - this will give fraction of bouts with
    % intervals > each value
    TempIntervalsBeforeFirstMotifSyll = cumsum(TempIntervalsBeforeFirstMotifSyll, 2);
    TempIntervalsBeforeFirstMotifSyll = TempIntervalsBeforeFirstMotifSyll > 0;
    
    AllBirdFractionofBoutsWithIntervalsGreaterThanVal(i,:) = sum(TempIntervalsBeforeFirstMotifSyll)/size(TempIntervalsBeforeFirstMotifSyll, 1);
end
   
% Now to plot the individual data and group data
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [416 333 900 650]);
p = panel();
p.pack({1});
p(1).select();
hold on;
plot(fliplr(Edges), 100*AllBirdFractionofBoutsWithIntervalsGreaterThanVal', 'b')
errorbar(fliplr(Edges), mean(100*AllBirdFractionofBoutsWithIntervalsGreaterThanVal), std(100*AllBirdFractionofBoutsWithIntervalsGreaterThanVal)/sqrt(size(AllBirdFractionofBoutsWithIntervalsGreaterThanVal,1)), 'k', 'LineWidth', 1);
p.fontsize = 16;
xlabel('Gap between syllables (msec)');
ylabel('Fraction of bouts with intervals > specified value (%)');
title('Intervals before the first motif syllable of a bout');
p.marginleft = 25;
p.marginright = 10;
p.margintop = 10;
p.marginbottom = 20;
axis([Edges(1) Edges(end) 0 100]);
plot(500 * ones(1,2), [0 100], 'k--');
plot([Edges(1) Edges(end)], 10 * ones(1,2), 'k--');
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(OutputDir, 'Inter-syllableGaps_BeforeFirstMotifSyllable.png'), '-dpng', '-r300');
  
%% Now next thing is to see how many intervals are greater than different values after the first motif syllable
Edges = 0:50:InterBoutInterval; % In ms
for i = 1:length(BirdParameters),
    ValidSongBouts = find((BirdParameters(i).Bouts(:,7) == 1) & (BirdParameters(i).Bouts(:,8) > 0) & (BirdParameters(i).Bouts(:,9) > 1));
    for j = 1:length(ValidSongBouts),
        FirstMotifSyll = NaN;
        for k = 1:length(BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Labels),
            if (~isempty(find(BirdParameters(i).MotifLabels == BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Labels(k))))
                FirstMotifSyll = k;
                break;
            end
        end
        if (isempty(find(isnan(FirstMotifSyll))))
            IntervalsAfterFirstMotifSyll{i}{j} = BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Onsets(FirstMotifSyll+1:end) - BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Offsets(FirstMotifSyll:end-1);
            BinnedIntervalsAfterFirstMotifSyll{i}(j,:) = histc(BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Onsets(FirstMotifSyll+1:end) - BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Offsets(FirstMotifSyll:end-1), Edges);
        end
    end
end

%% Now to find the fraction of bouts that have atleast 1 interval > each of the values in edges - sort of a cumulative fraction
for i = 1:length(BinnedIntervalsAfterFirstMotifSyll),
    TempIntervalsAfterFirstMotifSyll = BinnedIntervalsAfterFirstMotifSyll{i};
    % Threshold it to 1 and 0 for anything that is > 1
    TempIntervalsAfterFirstMotifSyll = TempIntervalsAfterFirstMotifSyll > 0;
    % Flip lr to get fraction of intervals > a particular value
    TempIntervalsAfterFirstMotifSyll = fliplr(TempIntervalsAfterFirstMotifSyll);
    % Now do a cumlative sum along the row and threshold it and divide by
    % total number of bouts - this will give fraction of bouts with
    % intervals > each value
    TempIntervalsAfterFirstMotifSyll = cumsum(TempIntervalsAfterFirstMotifSyll, 2);
    TempIntervalsAfterFirstMotifSyll = TempIntervalsAfterFirstMotifSyll > 0;
    
    AllBirdFractionofBoutsWithIntervalsGreaterThanValAfter(i,:) = sum(TempIntervalsAfterFirstMotifSyll)/size(TempIntervalsAfterFirstMotifSyll, 1);
end
   
% Now to plot the individual data and group data
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [416 333 900 650]);
p = panel();
p.pack({1});
p(1).select();
hold on;
plot(fliplr(Edges), 100*AllBirdFractionofBoutsWithIntervalsGreaterThanValAfter', 'b')
errorbar(fliplr(Edges), mean(100*AllBirdFractionofBoutsWithIntervalsGreaterThanValAfter), std(100*AllBirdFractionofBoutsWithIntervalsGreaterThanValAfter)/sqrt(size(AllBirdFractionofBoutsWithIntervalsGreaterThanValAfter,1)), 'k', 'LineWidth', 1);
p.fontsize = 16;
xlabel('Gap between syllables (msec)');
ylabel('Fraction of bouts with intervals > specified value (%)');
title('Intervals after the first motif syllable of a bout');
p.marginleft = 25;
p.marginright = 10;
p.margintop = 10;
p.marginbottom = 20;
axis([Edges(1) Edges(end) 0 100]);
plot(500 * ones(1,2), [0 100], 'k--');
plot([Edges(1) Edges(end)], 10 * ones(1,2), 'k--');
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(OutputDir, 'Inter-syllableGaps_AfterFirstMotifSyllable.png'), '-dpng', '-r300');

%% Now to find fraction of bouts where there are motif syllables after an interval of a specified value after the first motif syllable
% In essence find the intervals after the first motif syllable and the last
% motif syllable in a bout. 

Edges = 0:50:InterBoutInterval; % In ms
for i = 1:length(BirdParameters),
    ValidSongBouts = find((BirdParameters(i).Bouts(:,7) == 1) & (BirdParameters(i).Bouts(:,8) > 0) & (BirdParameters(i).Bouts(:,9) > 1));
    for j = 1:length(ValidSongBouts),
        FirstMotifSyll = NaN;
        LastMotifSyll = NaN;
        for k = 1:length(BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Labels),
            if (~isempty(find(BirdParameters(i).MotifLabels == BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Labels(k))))
                FirstMotifSyll = k;
                break;
            end
        end
        
        if (isempty(find(isnan(FirstMotifSyll))))
            for k = length(BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Labels):-1:FirstMotifSyll,
                if (~isempty(find(BirdParameters(i).MotifLabels == BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Labels(k))))
                    LastMotifSyll = k;
                    break;
                end
            end
            IntervalsBetweenFirstAndLastMotifSyll{i}{j} = BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Onsets(FirstMotifSyll+1:LastMotifSyll) - BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Offsets(FirstMotifSyll:LastMotifSyll-1);
            BinnedIntervalsBetweenFirstAndLastMotifSyll{i}(j,:) = histc(BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Onsets(FirstMotifSyll+1:LastMotifSyll) - BirdParameters(i).BoutDetails{ValidSongBouts(j)}.Offsets(FirstMotifSyll:LastMotifSyll-1), Edges);
        end
    end
end

%% Now to find the fraction of bouts that have atleast 1 interval > each of the values in edges - sort of a cumulative fraction
for i = 1:length(BinnedIntervalsBetweenFirstAndLastMotifSyll),
    TempIntervalsBetweenFirstAndLastMotifSyll = BinnedIntervalsBetweenFirstAndLastMotifSyll{i};
    % Threshold it to 1 and 0 for anything that is > 1
    TempIntervalsBetweenFirstAndLastMotifSyll = TempIntervalsBetweenFirstAndLastMotifSyll > 0;
    % Flip lr to get fraction of intervals > a particular value
    TempIntervalsBetweenFirstAndLastMotifSyll = fliplr(TempIntervalsBetweenFirstAndLastMotifSyll);
    % Now do a cumlative sum along the row and threshold it and divide by
    % total number of bouts - this will give fraction of bouts with
    % intervals > each value
    TempIntervalsBetweenFirstAndLastMotifSyll = cumsum(TempIntervalsBetweenFirstAndLastMotifSyll, 2);
    TempIntervalsBetweenFirstAndLastMotifSyll = TempIntervalsBetweenFirstAndLastMotifSyll > 0;
    
    AllBirdFractionofBoutsWithIntervalsBetweenFirstAndLast(i,:) = sum(TempIntervalsBetweenFirstAndLastMotifSyll)/size(TempIntervalsBetweenFirstAndLastMotifSyll, 1);
end
   
% Now to plot the individual data and group data
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [416 333 900 650]);
p = panel();
p.pack({1});
p(1).select();
hold on;
plot(fliplr(Edges), 100*AllBirdFractionofBoutsWithIntervalsBetweenFirstAndLast', 'b')
errorbar(fliplr(Edges), mean(100*AllBirdFractionofBoutsWithIntervalsBetweenFirstAndLast), std(100*AllBirdFractionofBoutsWithIntervalsBetweenFirstAndLast)/sqrt(size(AllBirdFractionofBoutsWithIntervalsBetweenFirstAndLast,1)), 'k', 'LineWidth', 1);
p.fontsize = 16;
xlabel('Gap between syllables (msec)');
ylabel('Fraction of bouts with intervals > specified value (%)');
title('Intervals between first and last motif syllable of a bout');
p.marginleft = 25;
p.marginright = 10;
p.margintop = 10;
p.marginbottom = 20;
axis([Edges(1) Edges(end) 0 100]);
plot(500 * ones(1,2), [0 100], 'k--');
plot([Edges(1) Edges(end)], 10 * ones(1,2), 'k--');
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(OutputDir, 'Inter-syllableGaps_BetweenFirstAndLastMotifSyllable.png'), '-dpng', '-r300');

disp('Finished');