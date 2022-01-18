function [AllBirdMeanBoutPosition_MotifProb] = Prasanth_PlotBoutPos_MotifProb(BirdParameters, MinTrialNo)

% ========= Position in bout and the probability of an Motif at that position
% =========================================================================

% First thing to analyse is the position of Motifs within a bout. Here, we
% will consider all bouts and calculate the probability of an Motif occuring
% at each of the positions within a bout. To test for significance, we will
% consider same set of labels for each bout, but scramble them up within
% a bout and now calculate the probability of an Motif at each position. We
% will do this 1000 time and then take the middle 95% percent as the
% confidence interval - anything outside this will be significant.

for i = 1:length(BirdParameters),
    % Find all valid song bouts
    % Valid song bout is one with at least one song syllable and 2s before
    % and after the bout
    ValidBouts = find((BirdParameters(i).Bouts(:,7) == 1) & (BirdParameters(i).Bouts(:,8) > 0) & (BirdParameters(i).Bouts(:,9) > 1));
    disp([BirdParameters(i).BirdName, ': ', num2str(length(ValidBouts)), ' valid song bouts']);
    
    Index = 1;
    for j = ValidBouts(:)',
        BoutLabels{i}{Index} = BirdParameters(i).NoteInfo{BirdParameters(i).Bouts(j,3)}.labels(BirdParameters(i).Bouts(j,1):BirdParameters(i).Bouts(j,2));
        Index = Index + 1;
    end
    
    % Now to calculate Motif probability at each position
    % First make a matrix with # of rows equal to number of bouts and # of
    % columns equal to longest bout length
    % This matrix is filled with NaNs
    BoutPosition_MotiforNot{i} = ones(length(BoutLabels{i}), max(cellfun(@length, BoutLabels{i}))) * NaN;
    
    % Now for each syllable in each bout, check if it is an Motif or not and
    % then put in a 1 or a 0 respectively in that position.
    for j = 1:length(BoutLabels{i}),
        for k = 1:length(BoutLabels{i}{j}),
            if (~isempty(find(BirdParameters(i).MotifLabels == BoutLabels{i}{j}(k))))
                BoutPosition_MotiforNot{i}(j, k) = 1;
            else
                BoutPosition_MotiforNot{i}(j, k) = 0;
            end
        end
    end
    
    for j = 1:size(BoutPosition_MotiforNot{i}, 2),
        BoutPosition_MotifProb{i}(j) = length(find(BoutPosition_MotiforNot{i}(:,j) == 1))/length(find(~isnan(BoutPosition_MotiforNot{i}(:,j))));
    end
    
    % Now to scramble up the bouts and do the same thing again with the
    % scrambled bouts - do this a total of 1000 times
    
    for RandomBoutLabels = 1:1000,
        TempBoutPosition_MotiforNot = ones(length(BoutLabels{i}), max(cellfun(@length, BoutLabels{i}))) * NaN;
        for j = 1:length(BoutLabels{i}),
            TempBoutLabels = BoutLabels{i}{j}(randperm(length(BoutLabels{i}{j})));
            for k = 1:length(TempBoutLabels),
                if (~isempty(find(BirdParameters(i).MotifLabels == TempBoutLabels(k))))
                    TempBoutPosition_MotiforNot(j, k) = 1;
                else
                    TempBoutPosition_MotiforNot(j, k) = 0;
                end
            end
        end
        
        for j = 1:size(TempBoutPosition_MotiforNot, 2),
            RandBoutPosition_MotifProb{i}(RandomBoutLabels, j) = length(find(TempBoutPosition_MotiforNot(:,j) == 1))/length(find(~isnan(TempBoutPosition_MotiforNot(:,j))));
        end
    end
    
    ConfidenceIntervals{i}(1,:) = prctile(RandBoutPosition_MotifProb{i}, 2.5);
    ConfidenceIntervals{i}(2,:) = prctile(RandBoutPosition_MotifProb{i}, 97.5);
    
    close all;
    figure;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [443 548 1111 400]);
    plot(BoutPosition_MotifProb{i}, 'ko-');
    hold on;
    plot(ConfidenceIntervals{i}', 'k--');
    axis tight;
    Temp = axis;
    Temp = [0.5 (size(BoutPosition_MotiforNot{i},2)+0.5) 0 1];
    xlabel('Syllable position within the bout');
    ylabel('Probability of syllable being an Motif');
    title([BirdParameters(i).BirdName, ': n=', num2str(length(ValidBouts)), ' song bouts']);
    set(gcf, 'PaperPositionMode', 'auto');
    print([BirdParameters(i).BirdName, '.BoutPos.vs.MotifProb.', num2str(BirdParameters(i).Interboutinterval), '.png'], '-dpng', '-r300');
end

% Now put all birds together, if the position is not significant is 0, if
% it is greater than chance, then it is 1 and if it is lesser than chance
% it is 0

AllBirdBoutPosition_MotifProb = ones(length(BoutPosition_MotifProb), max(cellfun(@length, BoutPosition_MotifProb))) * NaN;
AllBirdBoutPosition_MotifProb_RawValues = ones(length(BoutPosition_MotifProb), max(cellfun(@length, BoutPosition_MotifProb))) * NaN;
for i = 1:length(BoutPosition_MotifProb),
    AllBirdBoutPosition_MotifProb_RawValues(i,1:length(BoutPosition_MotifProb{i})) = BoutPosition_MotifProb{i};
    for j = 1:length(BoutPosition_MotifProb{i}),
        if (BoutPosition_MotifProb{i}(j) > ConfidenceIntervals{i}(2,j))
            AllBirdBoutPosition_MotifProb(i,j) = 1;
        else
            if (BoutPosition_MotifProb{i}(j) < ConfidenceIntervals{i}(1,j))
                AllBirdBoutPosition_MotifProb(i,j) = -1;
            else
                AllBirdBoutPosition_MotifProb(i,j) = 0;
            end
        end
    end
end

close all;
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [163 530 1000 400]);
hold on;
for i = 1:size(AllBirdBoutPosition_MotifProb,1),
    imagesc(1:1:length(BoutPosition_MotifProb{i}), i, AllBirdBoutPosition_MotifProb(i,1:length(BoutPosition_MotifProb{i})));
end
colormap([0 0 1; 0.9 0.9 0.5; 1 0 0]);
colorbar('Ticks', [-0.667 0 0.667], 'TickLabels', {'< chance', '= chance', '> chance'});
axis tight;
xlabel('Position within the bout');
ylabel('Bird #');
title(['Position in the bout vs. Motif significance (n=', num2str(size(AllBirdBoutPosition_MotifProb,1)), ') birds']);
set(gcf, 'PaperPositionMode', 'auto');
print(['AllBirds.BoutPos.vs.MotifProb.', num2str(BirdParameters(1).Interboutinterval), '.png'], '-dpng', '-r300');

for i = 1:size(AllBirdBoutPosition_MotifProb,2),
    AllBirdMeanBoutPosition_MotifProb(i,:) = [nanmean(AllBirdBoutPosition_MotifProb_RawValues(:,i)) nanstd(AllBirdBoutPosition_MotifProb_RawValues(:,i))/sqrt(length(find(~isnan(AllBirdBoutPosition_MotifProb_RawValues(:,i)))))];
end


close all;
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [163 530 1000 400]);
hold on;
errorbar(AllBirdMeanBoutPosition_MotifProb(:,1), AllBirdMeanBoutPosition_MotifProb(:,2), 'ko-', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'LineWidth', 2);
% plot(AllBirdBoutPosition_MotifProb_RawValues', 'k', 'MarkerSize', 4, 'LineWidth', 0.25, 'Color', [0.9 0.9 0.9]);
axis tight;
Temp = axis;
Temp = [0.5 35 0 1];
axis(Temp);
xlabel('Position within the bout');
ylabel('Probability of syllable being an Motif');
title(['Probability of Motifs in a bout (n=', num2str(size(AllBirdBoutPosition_MotifProb,1)), ') birds']);
set(gcf, 'PaperPositionMode', 'auto');
print(['AllBirds.BoutPos.vs.MeanMotifProb.', num2str(BirdParameters(1).Interboutinterval), '.png'], '-dpng', '-r300');