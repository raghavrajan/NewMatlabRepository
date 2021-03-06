function [MinDistances] = Prasanth_IndividualFeaturePCA(BirdParameters, MinTrialNo)

% ================= PCA Analysis ==========================================
% We want to do a PCA on all features of syllables and then compare the
% similarity of INs to all other categories of syllables - motif syllables,
% short calls and long calls. We will do this in two ways
% 1. Consider only vocalisations that are part of song bouts - for this we
% will include only birds that have enough short calls and long calls along
% with INs and motif syllables in their song bouts. Here we will not
% include only FF in the list of features for PCA
% 2. Consider short calls and long calls from non-song bouts also and this
% analysis, we will exclude both FF and amplitude from list of features for
% PCA, since amplitude varies as a function of the bird's position from the
% microphone.

Colours = 'rgbcmk';
Symbols = 'o+d^s<>.';

for i = 1:length(BirdParameters),
    % Find all valid song bouts
    % Valid song bout is one with at least one song syllable and 2s before
    % and after the bout
    ValidBouts = find((BirdParameters(i).Bouts(:,8) > 0) & (BirdParameters(i).Bouts(:,9) > 1));
    disp([BirdParameters(i).BirdName, ': ', num2str(length(ValidBouts)), ' valid song bouts']);
    
    % Now I have to find all the syllables that belong to these bouts and
    % keep only these in a separate array called ValidSyllables
    ValidBoutSyllableIndices{i} = [];
    ValidBoutSyllableBoutIndices{i} = [];
    ValidBoutSyllableCategories{i} = []; % whether it is a motif syllable, an IN, a short call or a long call - vocalization categories
    ValidBoutSyllableBoutCategories{i} = []; % whether it is part of a song bout or not - bout categories
    for j = ValidBouts(:)',
        Indices = find(BirdParameters(i).SyllableListBoutNum == j);
        ValidBoutSyllableIndices{i} = [ValidBoutSyllableIndices{i}; Indices(:)];
        ValidBoutSyllableBoutIndices{i} = [ValidBoutSyllableBoutIndices{i}; ones(size(Indices(:)))*j];
        ValidBoutSyllableBoutCategories{i} = [ValidBoutSyllableBoutCategories{i}; ones(size(Indices(:)))*BirdParameters(i).Bouts(j,7)];
        TempBoutLabels = char(BirdParameters(i).SyllableData(Indices,1));
        for k = 1:length(TempBoutLabels),
            if (~isempty(find(BirdParameters(i).MotifLabels == TempBoutLabels(k))))
                ValidBoutSyllableCategories{i}(end+1) = 0;
            else
                if (~isempty(find(BirdParameters(i).INLabels == TempBoutLabels(k))))
                    ValidBoutSyllableCategories{i}(end+1) = 1;
                else
                    if (~isempty(find(BirdParameters(i).LongcallLabels == TempBoutLabels(k))))
                        ValidBoutSyllableCategories{i}(end+1) = 2;
                    else
                        if (~isempty(find(BirdParameters(i).ShortcallLabels == TempBoutLabels(k))))
                            ValidBoutSyllableCategories{i}(end+1) = 3;
                        else
                            ValidBoutSyllableCategories{i}(end+1) = 4;
                        end
                    end
                end
            end
        end
    end
        
    % Now go through and check if all syllables have atleast a certain
    % number of renditions specified by the MinTrialNo argument.
    UniqueSylls = unique(char(BirdParameters(i).SyllableData(ValidBoutSyllableIndices{i}, 1)));
    SyllsToBeRemoved = [];
    for j = 1:length(UniqueSylls),
        Indices = (find(char(BirdParameters(i).SyllableData(ValidBoutSyllableIndices{i},1)) == UniqueSylls(j)));
        if (length(Indices) < MinTrialNo)
            SyllsToBeRemoved = [SyllsToBeRemoved; Indices(:)];
        end
    end
    if (~isempty(SyllsToBeRemoved))
        ValidBoutSyllableIndices{i}(SyllsToBeRemoved) = [];
        ValidBoutSyllableBoutIndices{i}(SyllsToBeRemoved) = [];
        ValidBoutSyllableCategories{i}(SyllsToBeRemoved) = []; 
        ValidBoutSyllableBoutCategories{i}(SyllsToBeRemoved) = []; 
    end
    
    % Now for part 1 of analysis, find all song bouts and check if all
    % categories of syllables are there
    SongBouts = (find(ValidBoutSyllableBoutCategories{i} == 1));
    AllSyllCategoriesFlag = 1; % 0 if all syllable categories are present and 1 if not
    for j = 0:3,
        if (isempty(find(ValidBoutSyllableCategories{i}(SongBouts) == j)))
            AllSyllCategoriesFlag = 0;
            break;
        end
    end
    if (AllSyllCategoriesFlag == 1)
        SongBouts = ValidBoutSyllableIndices{i}(find(ValidBoutSyllableBoutCategories{i} == 1));
        PCAFeatures = BirdParameters(i).SAPFeatsMatrix(SongBouts, 1:8);
        % Remove NanRows
        [NanRows, NanCols] = find(isnan(PCAFeatures));
        PCAFeatures(unique(NanRows),:) = NaN;
        disp(['Removed ', num2str(length(unique(NanRows))), ' nan rows']);
        % Normalize each column by mean and std
        PCAFeatures = (PCAFeatures - repmat(nanmean(PCAFeatures), size(PCAFeatures, 1), 1))./repmat(nanstd(PCAFeatures), size(PCAFeatures, 1), 1);
        
        [PCA_Coeff, PCA_Score, PCA_Latent, PCA_Tsquared, PCA_Explained] = pca(PCAFeatures);
        
        % First plot explained variance
        close all;
        figure;
        set(gcf, 'Color', 'w');
        plot(PCA_Explained, 'ko-');
        xlabel('Principal component #');
        ylabel('Variance explained (%)');
        title(BirdParameters(i).BirdName);
        set(gcf, 'PaperPositionMode', 'auto');
        print([BirdParameters(i).BirdName, '.PCA_VarianceExplained.', num2str(BirdParameters(1).Interboutinterval), '.png'], '-dpng', '-r300');


        % Now to calculate distance of INs from other syllables
        close all;
        figure;
        set(gcf, 'Color', 'w');
        set(gcf, 'Position', [680 288 850 650]);
        hold on;
        for j = 1:length(BirdParameters(i).INLabels),
            IN_Indices = find((char(BirdParameters(i).SyllableData(SongBouts,1)) == BirdParameters(i).INLabels(j)) & (~isnan(PCAFeatures(:,1))));
            EucDistance{i}{j} = [];
            PlotLegend = [];
            % first for motif syllables
            for k = 1:length(BirdParameters(i).MotifLabels),
                Syll_Indices = find((char(BirdParameters(i).SyllableData(SongBouts,1)) == BirdParameters(i).MotifLabels(k)) & (~isnan(PCAFeatures(:,1))));
                Distances = pdist2(PCA_Score(IN_Indices,1:4), PCA_Score(Syll_Indices,1:4));
                EucDistance{i}{j}(end+1,:) = [mean(Distances(:)) std(Distances(:)) 0 double(BirdParameters(i).MotifLabels(k))];
                if ((j == 1) && ~isempty(Syll_Indices))
                    plot(PCA_Score(Syll_Indices,1), PCA_Score(Syll_Indices, 2), [Colours(1), Symbols(mod(k, length(Symbols)) + 1)]);
                    PlotLegend{end+1} = [BirdParameters(i).MotifLabels(k), ' - motif'];
                end
            end
            
            % For INs
            for k = 1:length(BirdParameters(i).INLabels),
                Syll_Indices = find((char(BirdParameters(i).SyllableData(SongBouts,1)) == BirdParameters(i).INLabels(k)) & (~isnan(PCAFeatures(:,1))));
                Distances = pdist2(PCA_Score(IN_Indices,1:4), PCA_Score(Syll_Indices,1:4));
                EucDistance{i}{j}(end+1,:) = [mean(Distances(:)) std(Distances(:)) 1 double(BirdParameters(i).INLabels(k))];
                if ((j == 1) && ~isempty(Syll_Indices))
                    plot(PCA_Score(Syll_Indices,1), PCA_Score(Syll_Indices, 2), [Colours(2), Symbols(mod(k, length(Symbols)) + 1)]);
                    PlotLegend{end+1} = [BirdParameters(i).INLabels(k), ' - INs'];
                end
            end
            
            % For Long calls
            for k = 1:length(BirdParameters(i).LongcallLabels),
                Syll_Indices = find((char(BirdParameters(i).SyllableData(SongBouts,1)) == BirdParameters(i).LongcallLabels(k)) & (~isnan(PCAFeatures(:,1))));
                Distances = pdist2(PCA_Score(IN_Indices,1:4), PCA_Score(Syll_Indices,1:4));
                EucDistance{i}{j}(end+1,:) = [mean(Distances(:)) std(Distances(:)) 2 double(BirdParameters(i).LongcallLabels(k))];
                if ((j == 1) && ~isempty(Syll_Indices))
                    plot(PCA_Score(Syll_Indices,1), PCA_Score(Syll_Indices, 2), [Colours(3), Symbols(mod(k, length(Symbols)) + 1)]);
                    PlotLegend{end+1} = [BirdParameters(i).LongcallLabels(k), ' - long call'];
                end
            end
            
            % For short calls
            for k = 1:length(BirdParameters(i).ShortcallLabels),
                Syll_Indices = find((char(BirdParameters(i).SyllableData(SongBouts,1)) == BirdParameters(i).ShortcallLabels(k)) & (~isnan(PCAFeatures(:,1))));
                Distances = pdist2(PCA_Score(IN_Indices,1:4), PCA_Score(Syll_Indices,1:4));
                EucDistance{i}{j}(end+1,:) = [mean(Distances(:)) std(Distances(:)) 3 double(BirdParameters(i).ShortcallLabels(k))];
                if ((j == 1) && ~isempty(Syll_Indices))
                    plot(PCA_Score(Syll_Indices,1), PCA_Score(Syll_Indices, 2), [Colours(4), Symbols(mod(k, length(Symbols)) + 1)]);
                    PlotLegend{end+1} = [BirdParameters(i).ShortcallLabels(k), ' - short call'];
                end
            end
            if (j== 1)
                legend(PlotLegend);
                xlabel('PC 1');
                ylabel('PC 2');
                title(BirdParameters(i).BirdName);
            end
            
            % Now for each of the intro notes, we want to know what is the
            % shortest distance from this intro note to a syllable in each
            % of the three classes
            
            if ~(length(find(isnan(EucDistance{i}{j}(:,1)))) == size(EucDistance{i}{j},1))
                for k = 0:3,
                    if (k ~= 1)
                        Indices = find(EucDistance{i}{j}(:,3) == k);
                        [MinDist, MinDistIndex] = min(EucDistance{i}{j}(Indices,1));
                        MinEucDistance{i}{j}(k+1,:) = [EucDistance{i}{j}(Indices(MinDistIndex),:) i];
                    else
                        Indices = find(EucDistance{i}{j}(:,3) == k);
                        SelfEucDistance{i}{j}(1,:) = [EucDistance{i}{j}(Indices(j),:) i];
                        if (length(Indices) > 1)
                            Indices(j) = [];
                            [MinDist, MinDistIndex] = min(EucDistance{i}{j}(Indices,1));
                            MinEucDistance{i}{j}(k+1,:) = [EucDistance{i}{j}(Indices(MinDistIndex),:) i];
                        else
                            MinEucDistance{i}{j}(k+1,:) = ones(1,5)*NaN;
                        end
                    end
                end
            end
        end
        axis tight;
        Temp = axis;
        Temp = [1.02*Temp(1) 1.15*Temp(2) 1.02*Temp(3) 1.15*Temp(4)];
        axis(Temp);
        set(gcf, 'PaperPositionMode', 'auto');
        print([BirdParameters(i).BirdName, '.PC_Scores.', num2str(BirdParameters(1).Interboutinterval), '.png'], '-dpng', '-r300');
    end
end

% Now to plot minimum distance of INs to syllables of all categories
AllBirdMinEucDistance = [];
AllBirdSelfEucDistance = [];

for i = 1:length(MinEucDistance),
    for j = 1:length(MinEucDistance{i}),
        AllBirdMinEucDistance(end+1,:) = [MinEucDistance{i}{j}(:,1); MinEucDistance{i}{j}(1,end)];
        AllBirdSelfEucDistance(end+1,:) = [SelfEucDistance{i}{j}(1) SelfEucDistance{i}{j}(1,end)];
    end
end

MinDistances = [AllBirdMinEucDistance(:,1:end-1) AllBirdSelfEucDistance(:,1)];
% Plot the minimum distance of each IN from each of the syllable categories
figure;
hold on;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [234 547 1200 400]);
for i = 1:size(MinDistances,2),
    MinEucDistanceBar(i) = bar(i, nanmean(MinDistances(:,i)));
    set(MinEucDistanceBar(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 2);
    NumNotNanValues(i) = length(find(~isnan(MinDistances(:,i))));
end
plot(MinDistances', 'ko-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1.5);
errorbar(nanmean(MinDistances), nanstd(MinDistances)./sqrt(NumNotNanValues), 'ko-', 'LineWidth', 2, 'MarkerSize', 8);
set(gca, 'XTick', [1 2 3 4 5], 'XTickLabel', {'Motif syllables', 'Other INs', 'Long calls', 'Short calls', 'With same IN'});
ylabel('Mean distance');
title(['Distance between each IN and syllables of other categories in song bouts (n=', num2str(length(unique(AllBirdMinEucDistance(:,end)))), ' birds)']);
axis tight;
Temp = axis;
Temp = [0.25 5.75 0 1.02*Temp(4)];
axis(Temp);
set(gcf, 'PaperPositionMode', 'auto');
print(['AllBirds.INDistance_ToOtherCategories.', num2str(BirdParameters(1).Interboutinterval), '.png'], '-dpng', '-r300');

% Now to find the minimum distance for each IN and see what type of
% syllable it is. Excluding INs
for i = 1:size(MinDistances,1),
    [MinimumDistance, MinimumDistanceIndex(i)] = min(MinDistances(i,[1 3 4]));
end

MotifMatches = length(find(MinimumDistanceIndex == 1));
LongCallMatches = length(find(MinimumDistanceIndex == 2));
ShortCallMatches = length(find(MinimumDistanceIndex == 3));

figure; % All pie chart plotting copied from reference page for pie
set(gcf, 'Color', 'w');
set(gcf, 'Position', [111 420 1000 500]);

PieChartHandle = pie([MotifMatches LongCallMatches ShortCallMatches], [1 1 1]);
PieChartText = findobj(PieChartHandle, 'Type', 'text');
PercentValues = get(PieChartText, 'String');

PieChartLabelString = {['Motif syllables: ', num2str(100*(MotifMatches)/length(MinimumDistanceIndex)), '%']; ['     Long calls:', num2str(100*(LongCallMatches)/length(MinimumDistanceIndex)), '%']; ['    Short calls:', num2str((100*ShortCallMatches)/length(MinimumDistanceIndex)), '%']};
% CombinedPieChartLabelString = strcat(PieChartLabelString, PercentValues);
CombinedPieChartLabelString = PieChartLabelString;
OldExtents_cell = get(PieChartText, 'Extent');
OldExtents = cell2mat(OldExtents_cell);

PieChartText(1).String = CombinedPieChartLabelString(1);
PieChartText(2).String = CombinedPieChartLabelString(2);
PieChartText(3).String = CombinedPieChartLabelString(3);

NewExtents_cell = get(PieChartText, 'Extent');
NewExtents = cell2mat(NewExtents_cell);
WidthChange = NewExtents(:,3) - OldExtents(:,3);

SignValues = sign(OldExtents(:,1));
Offset = SignValues.*(WidthChange/2);

TextPositions_Cell = get(PieChartText, {'Position'});
TextPositions = cell2mat(TextPositions_Cell);
TextPositions(:,1) = TextPositions(:,1) + Offset;

PieChartText(1).Position = TextPositions(1,:);
PieChartText(2).Position = TextPositions(2,:);
PieChartText(3).Position = TextPositions(3,:);

title(['Nearest neighbour for each IN (excluding INs): n=', num2str(length(MinimumDistanceIndex)), ' INs from ', num2str(length(unique(AllBirdMinEucDistance(:,end)))), ' birds']);
set(gcf, 'PaperPositionMode', 'auto');
print(['AllBirds.INNearestNeighbour.', num2str(BirdParameters(1).Interboutinterval), '.png'], '-dpng', '-r300');

for i = 1:size(MinDistances,1),
    [MinimumDistance, MinimumDistanceIndexIncludingINs(i)] = min(MinDistances(i,1:4));
end
disp(MinimumDistanceIndexIncludingINs);
disp('Finished');
