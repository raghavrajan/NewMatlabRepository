function [BoutPosition_SyllDistribution] = Prasanth_PlotBoutPosFitExp(BirdParameters, MinTrialNo)

% ========= Position in bout and the syllable distribution ================

% I want to analyse the distribution of syllables at different positions
% within a bout. I will fit an exponential to the distibution. In addition
% I will correlate this distribution with the distribution over all
% positions.

for i = 1:length(BirdParameters),
    % Find all valid song bouts
    % Valid song bout is one with at least one song syllable and 2s before
    % and after the bout
    ValidBouts = find((BirdParameters(i).Bouts(:,7) == 1) & (BirdParameters(i).Bouts(:,8) > 0) & (BirdParameters(i).Bouts(:,9) > 1));
    disp([BirdParameters(i).BirdName, ': ', num2str(length(ValidBouts)), ' valid song bouts']);
    
    Index = 1;
    for j = ValidBouts(:)',
        BoutLabels{i}{Index} = BirdParameters(i).NoteInfo{BirdParameters(i).Bouts(j,3)}.labels(BirdParameters(i).Bouts(j,1):BirdParameters(i).Bouts(j,2));
        BoutOnsets{i}{Index} = BirdParameters(i).NoteInfo{BirdParameters(i).Bouts(j,3)}.onsets(BirdParameters(i).Bouts(j,1):BirdParameters(i).Bouts(j,2));
        BoutOffsets{i}{Index} = BirdParameters(i).NoteInfo{BirdParameters(i).Bouts(j,3)}.offsets(BirdParameters(i).Bouts(j,1):BirdParameters(i).Bouts(j,2));
        BoutSyllDurations{i}{Index} = BoutOffsets{i}{Index} - BoutOnsets{i}{Index};
        Index = Index + 1;
    end

    % Next to put all durations together in one big matrix with NaN's for
    % positions where there are no syllables in that bout. Each row will be
    % a bout and each column will be the duration of syllables at that
    % position in the bout
    BoutPosition_SyllDurations{i} = ones(length(BoutSyllDurations{i}), max(cellfun(@length, BoutSyllDurations{i}))) * NaN;
    
    for j = 1:length(BoutSyllDurations{i}),
        BoutPosition_SyllDurations{i}(j, 1:length(BoutSyllDurations{i}{j})) = BoutSyllDurations{i}{j};
    end
    
    % Now to find the number of bouts with syllables at each position across bouts
    NumNotNanValues = [];
    for j = 1:size(BoutPosition_SyllDurations{i},2),
        NumNotNanValues(j) = length(find(~isnan(BoutPosition_SyllDurations{i}(:,j))));
    end
    
    % Next to bin it and plot it
    Edges = 0:10:500;
    BoutPosition_SyllDistribution{i} = [];
    for j = 1:size(BoutPosition_SyllDurations{i},2),
        Indices = find(~isnan(BoutPosition_SyllDurations{i}(:,j)));
        if (length(Indices) >= MinTrialNo)
            BoutPosition_SyllDistribution{i}(j,:) = histc(BoutPosition_SyllDurations{i}(Indices,j), Edges)/length(Indices);
        else
            break;
        end
    end
    
    if (isempty(BoutPosition_SyllDistribution{i}))
        continue;
    end
    
    Indices = find(NumNotNanValues >= MinTrialNo);
    AllSylls = BoutPosition_SyllDurations{i}(:,Indices);
    AllSylls = AllSylls(find(~isnan(AllSylls)));
    if (~isempty(AllSylls))
        BoutPosition_SyllDistribution{i}(end+1,:) = histc(AllSylls(:), Edges)/length(AllSylls(:)); % overall distribution of all syllables irrespective of position
    else
        BoutPosition_SyllDistribution{i}(end+1,:) = ones(size(Edges))*NaN;
    end
    % Now to plot mean syllable duration
    Indices = find(NumNotNanValues >= MinTrialNo);
    MeanDurationFigure(i) = figure;
    hold on;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [282 400 1400 350]);
    errorbar(nanmean(BoutPosition_SyllDurations{i}(:,Indices)), nanstd(BoutPosition_SyllDurations{i}(:,Indices))./sqrt(NumNotNanValues(Indices)), 'bo-', 'LineWidth', 1.5);
    xlabel('Syllable position within the bout');
    ylabel('Mean syllable duration (msec)');
    title(BirdParameters(i).BirdName);
    
    % For testing significance, we will randomly change durations within
    % each bout and then repeat the procedure - will be done 10000 times
    % and we can use a p-value of 0.05
    
    Indices = find(NumNotNanValues >= MinTrialNo);
    TempBoutPos_SyllDurations = BoutPosition_SyllDurations{i}(:,Indices);
    clear RandomBoutPos_SyllDurations RandomMeanBoutPos_SyllDurations RandomBoutPosition_SyllDistribution Random_SyllDurHistCorr;
    
    for Reps = 1:10000,
        RandomBoutPos_SyllDurations = TempBoutPos_SyllDurations;
        for j = 1:size(TempBoutPos_SyllDurations, 1);
            Row_NotNaNValues = find(~isnan(TempBoutPos_SyllDurations(j,:)));
            RandomBoutPos_SyllDurations(j,Row_NotNaNValues) = RandomBoutPos_SyllDurations(j,randperm(length(Row_NotNaNValues)));
        end
        RandomMeanBoutPos_SyllDurations(Reps,:) = nanmean(RandomBoutPos_SyllDurations);

        % Now calculate correlations to the syllable duration histogram
        % based on this randomised syllable durations matrix
        for j = 1:size(RandomBoutPos_SyllDurations,2),
            Indices = find(~isnan(RandomBoutPos_SyllDurations(:,j)));
            RandomBoutPosition_SyllDistribution(j,:) = histc(RandomBoutPos_SyllDurations(Indices,j), Edges)/length(Indices);
            Random_SyllDurHistCorr(Reps, j) = (BoutPosition_SyllDistribution{i}(end,:) * RandomBoutPosition_SyllDistribution(j,:)')/(norm(BoutPosition_SyllDistribution{i}(end,:))*norm(RandomBoutPosition_SyllDistribution(j,:)));
        end
    end
    
    % Now plot the randomised mean 
    plot(mean(RandomMeanBoutPos_SyllDurations), 'ko-');
    % Now plot the 95% confidence intervals as dashed lines
    plot(prctile(RandomMeanBoutPos_SyllDurations, 97.5), 'k--');
    plot(prctile(RandomMeanBoutPos_SyllDurations, 2.5), 'k--');
    axis tight;
    Temp = axis;
    Temp = [0.5 size(RandomMeanBoutPos_SyllDurations, 2) 0.98*Temp(3) 1.05*Temp(4)];
    axis(Temp);
    
    % Now mark the ones that are significantly different from the random
    % calculated means
    Indices = find(NumNotNanValues >= MinTrialNo);
    SigIndices = zeros(size(Indices));
	for j = Indices(:)',
        if ((nanmean(BoutPosition_SyllDurations{i}(:,j)) > prctile(RandomMeanBoutPos_SyllDurations(:,j), 97.5)) || (nanmean(BoutPosition_SyllDurations{i}(:,j)) < prctile(RandomMeanBoutPos_SyllDurations(:,j), 2.5)))
            text(j, (Temp(4) * 1.03/1.05), '*', 'FontSize', 16);
            SigIndices(j) = 1;
        end
    end
    
    % Now to calculate and plot correlations between the histogram at each
    % position and the histogram for overall syllable durations across all
    % positions
    
    Indices = find(NumNotNanValues >= MinTrialNo);
    CorrSigIndices = zeros(size(Indices));
    for j = Indices(:)',
        SyllDurHistCorr{i}(j) = (BoutPosition_SyllDistribution{i}(end,:) * BoutPosition_SyllDistribution{i}(j,:)')/(norm(BoutPosition_SyllDistribution{i}(end,:))*norm(BoutPosition_SyllDistribution{i}(j,:)));
        if ((SyllDurHistCorr{i}(j) > prctile(Random_SyllDurHistCorr(:,j), 97.5)) || (SyllDurHistCorr{i}(j) < prctile(Random_SyllDurHistCorr(:,j), 2.5)))
            CorrSigIndices(j) = 1;
        end
    end
    CorrFigure(i) = figure;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [282 400 1400 350]);
    hold on;
    plot(SyllDurHistCorr{i}, 'bo-', 'LineWidth', 1.5);
    xlabel('Syllable position within the bout');
    ylabel('Correlation with overall syllable duration histogram');
    title(BirdParameters(i).BirdName);
    % Now plot the mean correlation for the random distribution
    plot(mean(Random_SyllDurHistCorr), 'ko-');
    % Now plot the 95% confidence intervals
    plot(prctile(Random_SyllDurHistCorr, 97.5), 'k--');
    plot(prctile(Random_SyllDurHistCorr, 2.5), 'k--');
    axis tight;
    Temp = axis;
    Temp = [0.5 size(Random_SyllDurHistCorr, 2) 0 1.05*Temp(4)];
    axis(Temp);
    for j = 1:size(Random_SyllDurHistCorr, 2),
        if (CorrSigIndices(j) == 1)
            text(j, Temp(4) * 1.03/1.05, '*', 'FontSize', 16);
        end
    end
    
    % Also, I will plot a matrix of correlations between each position with
    % the other
    Indices = find(NumNotNanValues >= MinTrialNo);
    for j = Indices(:)',
        for k = Indices(:)',
            AllSyllDurHistCorr{i}(j,k) = (BoutPosition_SyllDistribution{i}(j,:) * BoutPosition_SyllDistribution{i}(k,:)')/(norm(BoutPosition_SyllDistribution{i}(j,:))*norm(BoutPosition_SyllDistribution{i}(k,:)));
        end
    end
    figure;
    imagesc(AllSyllDurHistCorr{i});
    colorbar;
    xlabel('Syllable position within the bout');
    ylabel('Syllable position within the bout');
    
    RawSyllDistributionFigure(i) = figure;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [680 211 600 750]);
    hold on;
    for j = 1:size(BoutPosition_SyllDistribution{i}, 1),
        if (j ~= size(BoutPosition_SyllDistribution{i}, 1))
            if (CorrSigIndices(j) == 1)
                plot3(Edges, ones(size(Edges))*1000*(j+1), BoutPosition_SyllDistribution{i}(j,:), 'c', 'LineWidth', 3);
            else
                plot3(Edges, ones(size(Edges))*1000*(j+1), BoutPosition_SyllDistribution{i}(j,:), 'k', 'LineWidth', 1);
            end
        else
            plot3(Edges, ones(size(Edges))*0, BoutPosition_SyllDistribution{i}(j,:), 'r', 'LineWidth', 2);
        end
    end
    set(gca, 'YColor', 'w');
    set(gca, 'ZColor', 'w');
    view(0, 75);
    xlabel('Syll duration (ms)');
    title(BirdParameters(i).BirdName);
end

% Now plot the group data for syllable durations and for correlations

% First duration
AllBirdMeanBoutPos_SyllDurations = ones(length(BirdParameters), max(cellfun(@length, SyllDurHistCorr)))*NaN;

for i = 1:length(BoutPosition_SyllDurations),
    % Now to find the number of bouts with syllables at each position across bouts
    NumNotNanValues = [];
    for j = 1:size(BoutPosition_SyllDurations{i},2),
        NumNotNanValues(j) = length(find(~isnan(BoutPosition_SyllDurations{i}(:,j))));
    end
    Indices = find(NumNotNanValues >= MinTrialNo);
    AllBirdMeanBoutPos_SyllDurations(i,1:length(Indices)) = nanmean(BoutPosition_SyllDurations{i}(:,Indices));
end
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 297 1100 650]);
hold on;
plot(AllBirdMeanBoutPos_SyllDurations', 'ko-', 'LineWidth', 0.5);
for i = 1:size(AllBirdMeanBoutPos_SyllDurations, 2),
    MeanBoutPos_SyllDurBar(i) = bar(i, nanmean(AllBirdMeanBoutPos_SyllDurations(:,i)));
    set(MeanBoutPos_SyllDurBar(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
    Indices = find(~isnan(AllBirdMeanBoutPos_SyllDurations(:,i)));
    Mean_AllBirdBoutPosSyllDurations(i,:) = [nanmean(AllBirdMeanBoutPos_SyllDurations(:,i)) nanstd(AllBirdMeanBoutPos_SyllDurations(:,i))/sqrt(length(Indices))];
end
errorbar(Mean_AllBirdBoutPosSyllDurations(:,1), Mean_AllBirdBoutPosSyllDurations(:,2), 'ko-', 'LineWidth', 2, 'MarkerFaceColor', 'k');
axis tight;
Temp = axis;
Temp = [0.5 35.5 0 1.04*Temp(4)];
axis(Temp);
xlabel('Syllable position within the bout');
ylabel('Mean syllable duration (msec)');
title(['Mean syllable duration vs. position within the bout (n=', num2str(length(BirdParameters)), ' birds)']);
% Do the linear correlation for 1st 6 positions vs. syllable duration and
% then for the rest of the positions and syllable duration
PositionIndices = repmat(1:1:6, size(AllBirdMeanBoutPos_SyllDurations, 1), 1);
First6SyllDurations = AllBirdMeanBoutPos_SyllDurations(:,1:6);
[r, p] = corr(PositionIndices(:), First6SyllDurations(:), 'rows', 'complete');
disp(['Linear correlation between first six positions and mean syllable duration: r=', num2str(r), ';p=', num2str(p)]);

PositionIndices = repmat(6:1:15, size(AllBirdMeanBoutPos_SyllDurations, 1), 1);
RestSyllDurations = AllBirdMeanBoutPos_SyllDurations(:,6:15);
[r, p] = corr(PositionIndices(:), RestSyllDurations(:), 'rows', 'complete');
disp(['Linear correlation between 6th to 15th position and mean syllable duration: r=', num2str(r), ';p=', num2str(p)]);

% Next correlation
AllBirdMeanBoutPos_SyllCorr = ones(length(BirdParameters), max(cellfun(@length, SyllDurHistCorr)))*NaN;

for i = 1:length(SyllDurHistCorr),
    AllBirdMeanBoutPos_SyllCorr(i,1:length(SyllDurHistCorr{i})) = SyllDurHistCorr{i};
end
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 297 1100 650]);
hold on;
plot(AllBirdMeanBoutPos_SyllCorr', 'ko-', 'LineWidth', 0.5);
for i = 1:size(AllBirdMeanBoutPos_SyllCorr, 2),
    MeanBoutPos_SyllCorrBar(i) = bar(i, nanmean(AllBirdMeanBoutPos_SyllCorr(:,i)));
    set(MeanBoutPos_SyllCorrBar(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 1.5);
    Indices = find(~isnan(AllBirdMeanBoutPos_SyllCorr(:,i)));
    Mean_AllBirdBoutPosSyllCorr(i,:) = [nanmean(AllBirdMeanBoutPos_SyllCorr(:,i)) nanstd(AllBirdMeanBoutPos_SyllCorr(:,i))/sqrt(length(Indices))];
end
errorbar(Mean_AllBirdBoutPosSyllCorr(:,1), Mean_AllBirdBoutPosSyllCorr(:,2), 'ko-', 'LineWidth', 2, 'MarkerFaceColor', 'k');
axis tight;
Temp = axis;
Temp = [0.5 35.5 0 1.04*Temp(4)];
axis(Temp);
xlabel('Syllable position within the bout');
ylabel('Correlation with overall syllable duration histogram');
title(['n=', num2str(length(BirdParameters)), ' birds']);

PositionIndices = repmat(1:1:6, size(AllBirdMeanBoutPos_SyllCorr, 1), 1);
First6SyllCorr = AllBirdMeanBoutPos_SyllCorr(:,1:6);
[r, p] = corr(PositionIndices(:), First6SyllCorr(:), 'rows', 'complete');
disp(['Linear correlation between first six positions and correlation to overall syllable distribution histogram: r=', num2str(r), ';p=', num2str(p)]);

PositionIndices = repmat(6:1:15, size(AllBirdMeanBoutPos_SyllCorr, 1), 1);
RestSyllCorr = AllBirdMeanBoutPos_SyllCorr(:,6:15);
[r, p] = corr(PositionIndices(:), RestSyllCorr(:), 'rows', 'complete');
disp(['Linear correlation between 6th to 15th position and correlation to overall syllable distribution histogram: r=', num2str(r), ';p=', num2str(p)]);

disp('Finished');
%set(gcf, 'PaperPositionMode', 'auto');
%print(['AllBirds.BoutPos.vs.INProb.png'], '-dpng', '-r300');