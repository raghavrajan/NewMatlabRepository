function [Outputs] = Prasanth_AllBirdFeaturePCA(BirdParameters, MinTrialNo, varargin)

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
% This PCA analysis is done after combining syllables from all birds
% We will take all valid bouts and take all syllables from these bouts. Put
% them all together across birds and normalize the whole thing and then do
% PCA


Colours = 'rbkc';
Symbols = 'o+d^s<>.';
rng('default'); % To initialise random number generator

AllBirdSAPFeatures = [];
AllBirdSyllableCategories = [];
AllBirdSyllableLabels = [];
AllBirdSyllablePositionIndices = [];

for i = 1:length(BirdParameters),
    % Find all valid song bouts
    % Valid song bout is one with at least one song syllable and 2s before
    % and after the bout
    ValidBouts = find((BirdParameters(i).Bouts(:,8) > 0) & (BirdParameters(i).Bouts(:,9) > 1) & (BirdParameters(i).Bouts(:,7) == 1)); % only for song bouts with 2s before and after
    disp([BirdParameters(i).BirdName, ': ', num2str(length(ValidBouts)), ' valid song bouts']);
    
    % Now find all the valid syllables that are part of these bouts and get
    % their SAP features and their labels
    TempBoutLabels = [];
    TempSAPFeatures = [];
    TempSyllLabels = [];
    TempBoutPositionIndices = [];
    
    for j = ValidBouts(:)',
        Indices = find(BirdParameters(i).SyllableListBoutNum == j);
        TempSAPFeatures = [TempSAPFeatures; BirdParameters(i).SAPFeatsMatrix(Indices, 1:9)];
        TempSyllLabels = [TempSyllLabels; char(BirdParameters(i).SyllableData(Indices,1))];
        TempBoutLabels = [TempBoutLabels; char(BirdParameters(i).SyllableData(Indices,1))];
        TempBoutPositionIndices = [TempBoutPositionIndices; (1:1:length(Indices))'];
    end

    % I want to try and keep only 100 syllables of each type since I have a
    % lot of syllables
    if (nargin <= 2)
        TempUniqueSylls = unique(TempSyllLabels);
        SyllsToBeRetained = [];
        for j = 1:length(TempUniqueSylls),
            MatchingIndices = find(TempSyllLabels == TempUniqueSylls(j));
            [NanRows, NanCols] = find(isnan(TempSAPFeatures(MatchingIndices,:)));
            MatchingIndices(unique(NanRows)) = [];
            MatchingIndices = MatchingIndices(1:(min(length(MatchingIndices), 100)));
            SyllsToBeRetained = [SyllsToBeRetained; MatchingIndices];
        end
    else
        SyllsToBeRetained = 1:1:length(TempSyllLabels);
    end
    
    AllBirdSAPFeatures = [AllBirdSAPFeatures; TempSAPFeatures(SyllsToBeRetained,:)];
    AllBirdSyllableLabels = [AllBirdSyllableLabels; TempSyllLabels(SyllsToBeRetained)];
    TempBoutLabels = TempBoutLabels(SyllsToBeRetained);
    AllBirdSyllablePositionIndices = [AllBirdSyllablePositionIndices; TempBoutPositionIndices(SyllsToBeRetained)];
    
    for k = 1:length(TempBoutLabels),
        if (~isempty(find(BirdParameters(i).MotifLabels == TempBoutLabels(k))))
            AllBirdSyllableCategories(end+1,1) = 0;
        else
            if (~isempty(find(BirdParameters(i).INLabels == TempBoutLabels(k))))
                AllBirdSyllableCategories(end+1,1) = 1;
            else
                if (~isempty(find(BirdParameters(i).LongcallLabels == TempBoutLabels(k))))
                    AllBirdSyllableCategories(end+1,1) = 2;
                else
                    if (~isempty(find(BirdParameters(i).ShortcallLabels == TempBoutLabels(k))))
                        AllBirdSyllableCategories(end+1,1) = 3;
                    else
                        AllBirdSyllableCategories(end+1,1) = 4;
                    end
                end
            end
        end
    end
end

% Remove Nan values
[nanrows, nancols] = find(isnan(AllBirdSAPFeatures));
AllBirdSAPFeatures(unique(nanrows),:) = [];
AllBirdSyllableCategories(unique(nanrows),:) = [];
AllBirdSyllableLabels(unique(nanrows),:) = [];
AllBirdSyllablePositionIndices(unique(nanrows),:) = [];

disp(['Removed ', num2str(100 * length(unique(nanrows))/size(AllBirdSAPFeatures, 1)), '% of syllables with nan values']);

Index = 1;

ColourValues = [1 0 0; 0 0 1; 0 0 0; 0 1 1];

% rng('default'); % To start up the random number generator in the same way
for i = [2 3 4 6 7 9],
    figure;
    % subplot(3,2,Index);
    hold on;
    PlotLegend = [];
    for j = 0:3,
        Indices = find(AllBirdSyllableCategories == j);
        % Indices2 = Indices(randperm(min(5000, length(Indices))));
        Indices2 = Indices; % take all instead of a random subset
        scatter(AllBirdSAPFeatures(Indices2,1), AllBirdSAPFeatures(Indices2,i), ones(size(Indices2))*6, ones(length(Indices2), 1)*ColourValues(j+1,:), 'filled');
        switch (j)
            case 0
                PlotLegend{end+1} = 'Motif syllables';
            case 1
                PlotLegend{end+1} = 'INs';
            case 2
                PlotLegend{end+1} = 'Long calls';
            case 3
                PlotLegend{end+1} = 'Short calls';
        end
    end
    
    for j = 0:3,
        Indices = find(AllBirdSyllableCategories == j);
        if (~isempty(Indices))
            PlotConfidenceEllipse(AllBirdSAPFeatures(Indices, [1 i]), Colours(j+1), 1);
        end
    end
    
    %set(gca, 'Color', 'k');
    xlabel(BirdParameters(1).SAPFeat_FieldNames{1});
    ylabel(BirdParameters(1).SAPFeat_FieldNames{i});
    %legend(PlotLegend, 'TextColor', 'w');
    Index = Index + 1;
    axis tight;
end

% Do PCA on this data
if (nargin <= 2)
    [PCA_Coeff, PCA_Score, PCA_Latent, PCA_TSquared, PCA_Explained] = pca(zscore(AllBirdSAPFeatures));
    Outputs.PCA_Coeff = PCA_Coeff;
    Outputs.PCA_Explained = PCA_Explained;
else
    Inputs = varargin{1};
    PCA_Coeff = Inputs.PCA_Coeff;
    PCA_Explained = Inputs.PCA_Explained;
    AllBirdDataForEllipse = Inputs.AllBirdDataForEllipse;
    PCA_Score = zscore(AllBirdSAPFeatures)*PCA_Coeff;
    Outputs = Inputs;
end

% First plot explained variance
figure;
plot(PCA_Explained, 'ko-');
xlabel('Principal component #');
ylabel('Variance explained (%)');

% Next plot the first 2 PCs
figure;
set(gcf, 'Position', [218 472 1300 500]);
set(gcf, 'Color', 'w');
hold on;
PlotObjects = [];

for j = 0:3,
    Indices = find(AllBirdSyllableCategories == j);
%     if (nargin <= 2)
%         Indices2 = Indices(randperm(min(25000, length(Indices))));
%     else
%         Indices2 = Indices(randperm(min(75000, length(Indices))));
%     end
    Indices2 = Indices;
    
    subplot(1, 2, 1);
    hold on;
    % scatter(PCA_Score(Indices2,1), PCA_Score(Indices2,2), ones(size(Indices2))*6, ones(length(Indices2), 1)*ColourValues(j+1,:), 'filled');
    plot(PCA_Score(Indices2,1), PCA_Score(Indices2,2), [Colours(j+1), '.'], 'MarkerSize', 3 , 'MarkerFaceColor', ColourValues(j+1,:));
    % Now plot syllables from the first 3 positions onto this graph in
    % different colours
    if (nargin > 2)
        ExtraColours = 'cbk';
        for SyllPos = 1:2,
            Indices = find(AllBirdSyllablePositionIndices == SyllPos);
            plot(PCA_Score(Indices, 1), PCA_Score(Indices, 2), [ExtraColours(SyllPos), 'o'], 'MarkerSize', 2, 'MarkerFaceColor', ExtraColours(SyllPos));
        end
    end
    
    subplot(1, 2, 2);
    hold on;
    % PlotObjects(end+1) = scatter(PCA_Score(Indices2,1), PCA_Score(Indices2,3), ones(size(Indices2))*6, ones(length(Indices2), 1)*ColourValues(j+1,:), 'filled');
    if (nargin <= 2)
        PlotObjects(end+1) = plot(PCA_Score(Indices2,1), PCA_Score(Indices2,3), [Colours(j+1), '.'], 'MarkerSize', 3 , 'MarkerFaceColor', ColourValues(j+1,:));
    else
        plot(PCA_Score(Indices2,1), PCA_Score(Indices2,3), [Colours(j+1), '.'], 'MarkerSize', 3 , 'MarkerFaceColor', ColourValues(j+1,:));
    end
    if (nargin > 2)
        ExtraColours = 'cbk';
        for SyllPos = 1:2,
            Indices = find(AllBirdSyllablePositionIndices == SyllPos);
            plot(PCA_Score(Indices, 1), PCA_Score(Indices, 3), [ExtraColours(SyllPos), 'o'], 'MarkerSize', 2, 'MarkerFaceColor', ExtraColours(SyllPos));
        end
    end
    
    switch (j)
        case 0
            PlotLegend{end+1} = 'Motif syllables';
        case 1
            PlotLegend{end+1} = 'INs';
        case 2
            PlotLegend{end+1} = 'Long calls';
        case 3
            PlotLegend{end+1} = 'Short calls';
    end
end
for j = 0:3,
    if (nargin <= 2)
        Indices = find(AllBirdSyllableCategories == j);
        if (~isempty(Indices))
            subplot(1,2,1);
            PlotConfidenceEllipse(PCA_Score(Indices, [1 2]), Colours(j+1), 1);
            plot(mean(PCA_Score(Indices,1)), mean(PCA_Score(Indices,2)), [Colours(j+1), '+'], 'MarkerFaceColor', 'w', 'MarkerSize', 14, 'LineWidth', 3);
            
            subplot(1,2,2);
            PlotConfidenceEllipse(PCA_Score(Indices, [1 3]), Colours(j+1), 1);
            if (nargin <= 2)
                PlotObjects(end+1) = plot(mean(PCA_Score(Indices,1)), mean(PCA_Score(Indices,3)), [Colours(j+1), '+'], 'MarkerFaceColor', 'w', 'MarkerSize', 14, 'LineWidth', 3);
            else
                plot(mean(PCA_Score(Indices,1)), mean(PCA_Score(Indices,3)), [Colours(j+1), '+'], 'MarkerFaceColor', 'w', 'MarkerSize', 14, 'LineWidth', 3);
            end
        end
        Outputs.AllBirdDataForEllipse{j+1} = PCA_Score(Indices,:);
    else
        subplot(1, 2, 1);
        PlotConfidenceEllipse(AllBirdDataForEllipse{j+1}(:, [1 2]), Colours(j+1), 1);
        
        subplot(1, 2, 2);
        PlotConfidenceEllipse(AllBirdDataForEllipse{j+1}(:, [1 3]), Colours(j+1), 1);
    end
end
subplot(1, 2, 1);
axis tight;
Temp = axis;
Temp = 1.05*Temp;
axis(Temp);
xlabel('PC 1');
ylabel('PC 2');
if (nargin <= 2)
    legend(PlotObjects(1:4), PlotLegend(1:4));
end

subplot(1, 2, 2);
axis tight;
Temp = axis;
Temp = 1.05*Temp;
axis(Temp);
xlabel('PC 1');
ylabel('PC 3');
if (nargin <= 2)
    legend(PlotObjects(1:4), PlotLegend(1:4));
end

% Plot the co-effs
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [715 173 775 450]);
hold on;
for j = 1:4,
    plot(PCA_Coeff(:,j), [Colours(j), 'o-']);
    PCALegend{j} = ['PC #',num2str(j)];
end
plot([0 (size(PCA_Coeff, 2) + 0.5)], [0 0], 'k--');
set(gca, 'XTick', 1:1:size(PCA_Coeff, 2), 'XTickLabel', BirdParameters(1).SAPFeat_FieldNames(1:9), 'XTickLabelRotation', 45);
axis tight;
legend(PCALegend);
Temp = axis;
Temp = [0.5 (size(PCA_Coeff,2) + 3) 1.02*Temp(3) 1.02*Temp(4)];
axis(Temp);
ylabel('Co-efficients');

if (nargin <= 2)
    % This is normal bird data
    % Now use the first 4 PCs and calculate distances within and between
    % clusters
    FirstFourCategoryData = [];
    IndexArray = [];
    for j = 0:3,
        Indices = find(AllBirdSyllableCategories == j);
        NumSylls(j+1) = length(Indices);
        IndexArray = [IndexArray; ones(length(Indices),1)*j];
        FirstFourCategoryData = [FirstFourCategoryData; PCA_Score(Indices,:)];
        
        for k = 0:3,
            Indices2 = find(AllBirdSyllableCategories == k);
            Distances{j+1,k+1} = pdist2(mean(PCA_Score(Indices,1:4)), PCA_Score(Indices2,1:4));
            MeanDistances(j+1,k+1) = mean(Distances{j+1,k+1}(:));
            STDDistances(j+1,k+1) = std(Distances{j+1,k+1}(:));
            SEMDistances(j+1,k+1) = STDDistances(j+1,k+1)/sqrt(length(Distances{j+1,k+1}));
        end
    end

    % Now to check what is the minimum distance based on random assignment
    % of syllables to the different categories
    for Reps = 1:10000,
        RandomData = FirstFourCategoryData(randperm(sum(NumSylls)), :);
        for j = 0:3,
            Indices = find(IndexArray == j);
            for k = 0:3,
                Indices2 = find(IndexArray == k);
                RDistances = pdist2(mean(FirstFourCategoryData(Indices,1:4)), RandomData(Indices2,1:4));
                RandMeanDistances{j+1}{k+1}(Reps) = mean(RDistances(:));
                RandSTDDistances{j+1}{k+1}(Reps) = std(RDistances(:));
                RandSEMDistances{j+1}{k+1}(Reps) = RandSTDDistances{j+1}{k+1}(Reps)/sqrt(length(RDistances));
            end
        end
    end
    figure;
    set(gcf, 'Color', 'w');
    hold on;
    errorbar(MeanDistances(2,[2 4 1 3]), SEMDistances(2,[2 4 1 3]), 'ko-', 'MarkerFaceColor', 'k');
    Temp = cellfun(@prctile, RandMeanDistances{2}, {2.5 2.5 2.5 2.5});
    plot(Temp([2 4 1 3]), 'k--', 'Color', [0.7 0.7 0.7]);
    Temp = cellfun(@prctile, RandMeanDistances{2}, {97.5 97.5 97.5 97.5});
    plot(Temp([2 4 1 3]), 'k--', 'Color', [0.7 0.7 0.7]);
    set(gca, 'XTick', [1 2 3 4], 'XTickLabel', {'INs' 'Short calls' 'Motif syllables', 'Long calls'});
    ylabel('Mean distance from centroid of IN cluster');
    
    figure;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [680 430 550 550]);
    hold on;
    errorbar(MeanDistances([2 4 1 3], 2), SEMDistances([2 4 1 3], 2), 'ko-', 'MarkerFaceColor', 'k');
    Temp = [prctile(RandMeanDistances{1}{2}, 2.5) prctile(RandMeanDistances{2}{2}, 2.5) prctile(RandMeanDistances{3}{2}, 2.5) prctile(RandMeanDistances{4}{2}, 2.5)];
    plot(Temp([2 4 1 3]), 'k--', 'Color', [0.7 0.7 0.7]);
    Temp = [prctile(RandMeanDistances{1}{2}, 97.5) prctile(RandMeanDistances{2}{2}, 97.5) prctile(RandMeanDistances{3}{2}, 97.5) prctile(RandMeanDistances{4}{2}, 97.5)];
    plot(Temp([2 4 1 3]), 'k--', 'Color', [0.7 0.7 0.7]);
    set(gca, 'XTick', [1 2 3 4], 'XTickLabel', {'INs' 'Short calls' 'Motif syllables', 'Long calls'});
    ylabel('Mean distance of INs from centroid of syllable type cluster');
else
    % Calculate the distances of first and second position syllables for
    % lesion bird data from the mean of the syllable categories for normal
    % data
    for SyllPos = 1:10,
        Indices = find(AllBirdSyllablePositionIndices == SyllPos);
        for j = 0:3,
            TempDistances = pdist2(PCA_Score(Indices,1:4), mean(AllBirdDataForEllipse{j+1}(:,1:4)));
            MeanFirst10SyllDistances(SyllPos, j+1) = mean(TempDistances(:));
            SEMFirst10SyllDistances(SyllPos, j+1) = std(TempDistances(:))/sqrt(length(Indices));
            % As a measure of random distances, I will take a random subset of
            % syllables and measure their distance relative to these 4
            % categories
            for Reps = 1:10000,
                RandomIndices = randperm(size(PCA_Score,1));
                RandomDistances = pdist2(PCA_Score(RandomIndices(1:length(Indices)), 1:4), mean(AllBirdDataForEllipse{j+1}(:,1:4)));
                RandomFirst10SyllDistances{SyllPos}(Reps, j+1) = mean(RandomDistances(:));
            end
        end
    end
end
disp('Finished');
