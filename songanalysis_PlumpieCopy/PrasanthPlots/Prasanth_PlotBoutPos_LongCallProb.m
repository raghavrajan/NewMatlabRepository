function [AllBirdMeanBoutPosition_LongcallProb] = Prasanth_PlotBoutPos_LongCallProb(BirdParameters, MinTrialNo)

% ========= Position in bout and the probability of an Longcall at that position
% =========================================================================

% First thing to analyse is the position of Longcalls within a bout. Here, we
% will consider all bouts and calculate the probability of an Longcall occuring
% at each of the positions within a bout. To test for significance, we will
% consider same set of labels for each bout, but scramble them up within
% a bout and now calculate the probability of an Longcall at each position. We
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
    
    % Now to calculate Longcall probability at each position
    % First make a matrix with # of rows equal to number of bouts and # of
    % columns equal to longest bout length
    % This matrix is filled with NaNs
    BoutPosition_LongcallorNot{i} = ones(length(BoutLabels{i}), max(cellfun(@length, BoutLabels{i}))) * NaN;
    
    % Now for each syllable in each bout, check if it is an Longcall or not and
    % then put in a 1 or a 0 respectively in that position.
    for j = 1:length(BoutLabels{i}),
        for k = 1:length(BoutLabels{i}{j}),
            if (~isempty(find(BirdParameters(i).LongcallLabels == BoutLabels{i}{j}(k))))
                BoutPosition_LongcallorNot{i}(j, k) = 1;
            else
                BoutPosition_LongcallorNot{i}(j, k) = 0;
            end
        end
    end
    
    for j = 1:size(BoutPosition_LongcallorNot{i}, 2),
        BoutPosition_LongcallProb{i}(j) = length(find(BoutPosition_LongcallorNot{i}(:,j) == 1))/length(find(~isnan(BoutPosition_LongcallorNot{i}(:,j))));
    end
    
    % Now to scramble up the bouts and do the same thing again with the
    % scrambled bouts - do this a total of 1000 times
    
    for RandomBoutLabels = 1:1000,
        TempBoutPosition_LongcallorNot = ones(length(BoutLabels{i}), max(cellfun(@length, BoutLabels{i}))) * NaN;
        for j = 1:length(BoutLabels{i}),
            TempBoutLabels = BoutLabels{i}{j}(randperm(length(BoutLabels{i}{j})));
            for k = 1:length(TempBoutLabels),
                if (~isempty(find(BirdParameters(i).LongcallLabels == TempBoutLabels(k))))
                    TempBoutPosition_LongcallorNot(j, k) = 1;
                else
                    TempBoutPosition_LongcallorNot(j, k) = 0;
                end
            end
        end
        
        for j = 1:size(TempBoutPosition_LongcallorNot, 2),
            RandBoutPosition_LongcallProb{i}(RandomBoutLabels, j) = length(find(TempBoutPosition_LongcallorNot(:,j) == 1))/length(find(~isnan(TempBoutPosition_LongcallorNot(:,j))));
        end
    end
    
    ConfidenceIntervals{i}(1,:) = prctile(RandBoutPosition_LongcallProb{i}, 2.5);
    ConfidenceIntervals{i}(2,:) = prctile(RandBoutPosition_LongcallProb{i}, 97.5);
    
    close all;
    figure;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [443 548 1111 400]);
    plot(BoutPosition_LongcallProb{i}, 'ko-');
    hold on;
    plot(ConfidenceIntervals{i}', 'k--');
    axis tight;
    Temp = axis;
    Temp = [0.5 (size(BoutPosition_LongcallorNot{i},2)+0.5) 0 1];
    xlabel('Syllable position within the bout');
    ylabel('Probability of syllable being an Longcall');
    title([BirdParameters(i).BirdName, ': n=', num2str(length(ValidBouts)), ' song bouts']);
    set(gcf, 'PaperPositionMode', 'auto');
    print([BirdParameters(i).BirdName, '.BoutPos.vs.LongcallProb.', num2str(BirdParameters(i).Interboutinterval), '.png'], '-dpng', '-r300');
end

% Now put all birds together, if the position is not significant is 0, if
% it is greater than chance, then it is 1 and if it is lesser than chance
% it is 0

AllBirdBoutPosition_LongcallProb = ones(length(BoutPosition_LongcallProb), max(cellfun(@length, BoutPosition_LongcallProb))) * NaN;
AllBirdBoutPosition_LongcallProb_RawValues = ones(length(BoutPosition_LongcallProb), max(cellfun(@length, BoutPosition_LongcallProb))) * NaN;
for i = 1:length(BoutPosition_LongcallProb),
    AllBirdBoutPosition_LongcallProb_RawValues(i,1:length(BoutPosition_LongcallProb{i})) = BoutPosition_LongcallProb{i};
    for j = 1:length(BoutPosition_LongcallProb{i}),
        if (BoutPosition_LongcallProb{i}(j) > ConfidenceIntervals{i}(2,j))
            AllBirdBoutPosition_LongcallProb(i,j) = 1;
        else
            if (BoutPosition_LongcallProb{i}(j) < ConfidenceIntervals{i}(1,j))
                AllBirdBoutPosition_LongcallProb(i,j) = -1;
            else
                AllBirdBoutPosition_LongcallProb(i,j) = 0;
            end
        end
    end
end

close all;
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [163 530 1000 400]);
hold on;
for i = 1:size(AllBirdBoutPosition_LongcallProb,1),
    imagesc(1:1:length(BoutPosition_LongcallProb{i}), i, AllBirdBoutPosition_LongcallProb(i,1:length(BoutPosition_LongcallProb{i})));
end
colormap([0 0 1; 0.9 0.9 0.5; 1 0 0]);
colorbar('Ticks', [-0.667 0 0.667], 'TickLabels', {'< chance', '= chance', '> chance'});
axis tight;
xlabel('Position within the bout');
ylabel('Bird #');
title(['Position in the bout vs. Longcall significance (n=', num2str(size(AllBirdBoutPosition_LongcallProb,1)), ') birds']);
set(gcf, 'PaperPositionMode', 'auto');
print(['AllBirds.BoutPos.vs.LongcallProb.', num2str(BirdParameters(1).Interboutinterval), '.png'], '-dpng', '-r300');

for i = 1:size(AllBirdBoutPosition_LongcallProb,2),
    AllBirdMeanBoutPosition_LongcallProb(i,:) = [nanmean(AllBirdBoutPosition_LongcallProb_RawValues(:,i)) nanstd(AllBirdBoutPosition_LongcallProb_RawValues(:,i))/sqrt(length(find(~isnan(AllBirdBoutPosition_LongcallProb_RawValues(:,i)))))];
end


close all;
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [163 530 1000 400]);
hold on;
errorbar(AllBirdMeanBoutPosition_LongcallProb(:,1), AllBirdMeanBoutPosition_LongcallProb(:,2), 'ko-', 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'LineWidth', 2);
% plot(AllBirdBoutPosition_LongcallProb_RawValues', 'k', 'MarkerSize', 4, 'LineWidth', 0.25, 'Color', [0.9 0.9 0.9]);
axis tight;
Temp = axis;
Temp = [0.5 35 0 1];
axis(Temp);
xlabel('Position within the bout');
ylabel('Probability of syllable being an Short call');
title(['Probability of long calls within a bout (n=', num2str(size(AllBirdBoutPosition_LongcallProb,1)), ') birds']);
set(gcf, 'PaperPositionMode', 'auto');
print(['AllBirds.BoutPos.vs.MeanLongcallProb.', num2str(BirdParameters(1).Interboutinterval), '.png'], '-dpng', '-r300');