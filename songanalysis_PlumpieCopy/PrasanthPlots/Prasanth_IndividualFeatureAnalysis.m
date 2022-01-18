function [] = Prasanth_IndividualFeatureAnalysis(BirdParameters, MinTrialNo)

% ================= Feature analysis - comparision of features for INs and
% other syllables =========================================================

% We want to compare the features of INs and all other syllables and see if
% INs in general have some common characteristics. Also, previous studies
% suggest a difference between INs and the motif and we want to quantify
% this

% First, we'll take bouts with song in it and take INs and motif syllables
% from that. Then we'll take calls from non-song bouts and include them in
% the analysis. First we'll check for each SAP feature individually and do
% the statistics on it.

% Now, I already have a list of syllables and the bouts where they belong.
% Now, I need to identify valid bouts (i.e. bouts with 2s before and after)
% and then keep only syllables belonging to this bouts. 

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
    % Now if we have to make separate figures for the following:
    % 1) Figure with all SAP features for all syllables separately
    % depending on whether they are part of a song bout or a non-song
    % bout
    % 2) Figure will all SAP features based on the different syllable
    % categories

    % Now, first figure (1) above
    close all;
    figure;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [362 69 1200 900]);
    for j = 1:length(BirdParameters(i).SAPFeat_FieldNames),
        subplot(ceil(length(BirdParameters(i).SAPFeat_FieldNames)/2), 2, j);
        hold on;
        PlotXLabelString = [];
        Index = 1;
        for k = 1:length(BirdParameters(i).MotifLabels),
            Indices = ValidBoutSyllableIndices{i}(find(char(BirdParameters(i).SyllableData(ValidBoutSyllableIndices{i},1)) == BirdParameters(i).MotifLabels(k)));
            if (~isempty(Indices))
                errorbar(Index, nanmean(BirdParameters(i).SAPFeatsMatrix(Indices,j)), nanstd(BirdParameters(i).SAPFeatsMatrix(Indices,j)), 'ko');
                PlotXLabelString{Index} = BirdParameters(i).MotifLabels(k);
            end
            Index = Index + 1;
        end

        Index = Index + 1;
        for k = 1:length(BirdParameters(i).INLabels),
            % First for song bouts
            Indices = ValidBoutSyllableIndices{i}(find((char(BirdParameters(i).SyllableData(ValidBoutSyllableIndices{i},1)) == BirdParameters(i).INLabels(k)) & (ValidBoutSyllableBoutCategories{i} == 1)));
            if (~isempty(Indices))
                errorbar(Index, nanmean(BirdParameters(i).SAPFeatsMatrix(Indices,j)), nanstd(BirdParameters(i).SAPFeatsMatrix(Indices,j)), 'ko');
                PlotXLabelString{Index} = ['SB-',BirdParameters(i).INLabels(k)];
            end
            Index = Index + 1;

            % Next for non-song bouts
            Indices = ValidBoutSyllableIndices{i}(find((char(BirdParameters(i).SyllableData(ValidBoutSyllableIndices{i},1)) == BirdParameters(i).INLabels(k)) & (ValidBoutSyllableBoutCategories{i} == 0)));
            if (~isempty(Indices))
                errorbar(Index, nanmean(BirdParameters(i).SAPFeatsMatrix(Indices,j)), nanstd(BirdParameters(i).SAPFeatsMatrix(Indices,j)), 'ko');
                PlotXLabelString{Index} = ['NSB-',BirdParameters(i).INLabels(k)];
            end
            Index = Index + 1;                
        end

        Index = Index + 1;
        for k = 1:length(BirdParameters(i).LongcallLabels),
            % First for song bouts
            Indices = ValidBoutSyllableIndices{i}(find((char(BirdParameters(i).SyllableData(ValidBoutSyllableIndices{i},1)) == BirdParameters(i).LongcallLabels(k)) & (ValidBoutSyllableBoutCategories{i} == 1)));
            if (~isempty(Indices))
                errorbar(Index, nanmean(BirdParameters(i).SAPFeatsMatrix(Indices,j)), nanstd(BirdParameters(i).SAPFeatsMatrix(Indices,j)), 'ko');
                PlotXLabelString{Index} = ['SB-',BirdParameters(i).LongcallLabels(k)];
            end
            Index = Index + 1;

            % Next for non-song bouts
            Indices = ValidBoutSyllableIndices{i}(find((char(BirdParameters(i).SyllableData(ValidBoutSyllableIndices{i},1)) == BirdParameters(i).LongcallLabels(k)) & (ValidBoutSyllableBoutCategories{i} == 0)));
            if (~isempty(Indices))
                errorbar(Index, nanmean(BirdParameters(i).SAPFeatsMatrix(Indices,j)), nanstd(BirdParameters(i).SAPFeatsMatrix(Indices,j)), 'ko');
                PlotXLabelString{Index} = ['NSB-',BirdParameters(i).LongcallLabels(k)];
            end
            Index = Index + 1;                
        end

        Index = Index + 1;
        for k = 1:length(BirdParameters(i).ShortcallLabels),
            % First for song bouts
            Indices = ValidBoutSyllableIndices{i}(find((char(BirdParameters(i).SyllableData(ValidBoutSyllableIndices{i},1)) == BirdParameters(i).ShortcallLabels(k)) & (ValidBoutSyllableBoutCategories{i} == 1)));
            if (~isempty(Indices))
                errorbar(Index, nanmean(BirdParameters(i).SAPFeatsMatrix(Indices,j)), nanstd(BirdParameters(i).SAPFeatsMatrix(Indices,j)), 'ko');
                PlotXLabelString{Index} = ['SB-',BirdParameters(i).ShortcallLabels(k)];
            end
            Index = Index + 1;

            % Next for non-song bouts
            Indices = ValidBoutSyllableIndices{i}(find((char(BirdParameters(i).SyllableData(ValidBoutSyllableIndices{i},1)) == BirdParameters(i).ShortcallLabels(k)) & (ValidBoutSyllableBoutCategories{i} == 0)));
            if (~isempty(Indices))
                errorbar(Index, nanmean(BirdParameters(i).SAPFeatsMatrix(Indices,j)), nanstd(BirdParameters(i).SAPFeatsMatrix(Indices,j)), 'ko');
                PlotXLabelString{Index} = ['NSB-',BirdParameters(i).ShortcallLabels(k)];
            end
            Index = Index + 1;                
        end
        
        ylabel(BirdParameters(i).SAPFeat_FieldNames{j});
        if (j > 6)
            set(gca, 'XTick', 1:1:length(PlotXLabelString), 'XTickLabel', PlotXLabelString, 'XTickLabelRotation', 45);
        else
            set(gca, 'XTick', 1:1:length(PlotXLabelString), 'XTickLabel', []);
        end
        if (j < 3)
            title(BirdParameters(i).BirdName);
        end
    end
    set(gcf, 'PaperPositionMode', 'auto');
    print([BirdParameters(i).BirdName, '.SyllFeaturesForIndividualSyllables.png'], '-dpng', '-r300');
    
    % Now to plot by categories
    close all;
    figure;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [362 69 1200 900]);
    for j = 1:length(BirdParameters(i).SAPFeat_FieldNames),
        subplot(ceil(length(BirdParameters(i).SAPFeat_FieldNames)/2), 2, j);
        hold on;
        PlotXLabelString = [];
        Index = 1;
        for k = min(ValidBoutSyllableCategories{i}):max(ValidBoutSyllableCategories{i}),
            % First for song bouts
            Indices = ValidBoutSyllableIndices{i}(find((ValidBoutSyllableCategories{i}(:) == k) & (ValidBoutSyllableBoutCategories{i}(:) == 1)));
            if (~isempty(Indices))
                errorbar(Index, nanmean(BirdParameters(i).SAPFeatsMatrix(Indices,j)), nanstd(BirdParameters(i).SAPFeatsMatrix(Indices,j)), 'ko');
                PlotXLabelString{Index} = ['SB - Syll Type ',num2str(k)];
            end
            Index = Index + 1;
            
            % Next for non-song bouts
            Indices = ValidBoutSyllableIndices{i}(find((ValidBoutSyllableCategories{i}(:) == k) & (ValidBoutSyllableBoutCategories{i}(:) == 0)));
            if (~isempty(Indices))
                errorbar(Index, nanmean(BirdParameters(i).SAPFeatsMatrix(Indices,j)), nanstd(BirdParameters(i).SAPFeatsMatrix(Indices,j)), 'ko');
                PlotXLabelString{Index} = ['NSB - Syll Type ',num2str(k)];
            end
            Index = Index + 1;                
        end

        ylabel(BirdParameters(i).SAPFeat_FieldNames{j});
        if (j > 6)
            set(gca, 'XTick', 1:1:length(PlotXLabelString), 'XTickLabel', PlotXLabelString, 'XTickLabelRotation', 45);
        else
            set(gca, 'XTick', 1:1:length(PlotXLabelString), 'XTickLabel', []);
        end        
        if (j < 3)
            title(BirdParameters(i).BirdName);
        end
    end
    set(gcf, 'PaperPositionMode', 'auto');
    print([BirdParameters(i).BirdName, '.SyllFeaturesForSyllCategories.png'], '-dpng', '-r300');
    
    % Now to plot only for song bouts and also put significances - in this
    % case we will only do the comparisons with INs and each of the other
    % categories.
    close all;
    figure;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [362 69 1200 900]);
    for j = 1:length(BirdParameters(i).SAPFeat_FieldNames),
        subplot(ceil(length(BirdParameters(i).SAPFeat_FieldNames)/2), 2, j);
        hold on;
        PlotXLabelString = [];
        DataToCompare = [];
        for k = 0:1:3,
            % First for song bouts
            Indices = (find((ValidBoutSyllableCategories{i}(:) == k) & (ValidBoutSyllableBoutCategories{i}(:) == 1)));
            if (length(Indices) >= MinTrialNo)
                DataToCompare = [DataToCompare; [BirdParameters(i).SAPFeatsMatrix(ValidBoutSyllableIndices{i}(Indices),j) ones(length(Indices),1)*k]];
                errorbar(k, nanmean(BirdParameters(i).SAPFeatsMatrix(ValidBoutSyllableIndices{i}(Indices),j)), nanstd(BirdParameters(i).SAPFeatsMatrix(ValidBoutSyllableIndices{i}(Indices),j)), 'ko');
                switch k
                    case 0
                        PlotXLabelString{k+1} = 'Motif syllables';
                    case 1
                        PlotXLabelString{k+1} = 'INs';
                    case 2
                        PlotXLabelString{k+1} = 'Long calls';
                    case 3
                        PlotXLabelString{k+1} = 'Short calls';
                end
            end
        end

        % Now do a anova comparison of all these 4 categories
        [p, tabl, stats] = anova1(DataToCompare(:,1), DataToCompare(:,2), 'off');
        Comparisons = multcompare(stats, 'display', 'off');
        INComparisonIndices = find((Comparisons(:,1) == 2) | (Comparisons(:,2) == 2));
        SignificanceValues = [];
        
        for k = INComparisonIndices(:)',
            if (Comparisons(k,end) < 0.05)
                if (Comparisons(k,1) == 2)
                    Indices = (find((ValidBoutSyllableCategories{i}(:) == (Comparisons(k,2)-1)) & (ValidBoutSyllableBoutCategories{i}(:) == 1)));
                    text(Comparisons(k,2) - 1, 1.01*(nanmean(BirdParameters(i).SAPFeatsMatrix(ValidBoutSyllableIndices{i}(Indices),j)) + nanstd(BirdParameters(i).SAPFeatsMatrix(ValidBoutSyllableIndices{i}(Indices),j))), '*', 'FontSize', 14);
                else
                    Indices = (find((ValidBoutSyllableCategories{i}(:) == (Comparisons(k,1)-1)) & (ValidBoutSyllableBoutCategories{i}(:) == 1)));
                    text(Comparisons(k,1) - 1, 1.05*(nanmean(BirdParameters(i).SAPFeatsMatrix(ValidBoutSyllableIndices{i}(Indices),j)) + nanstd(BirdParameters(i).SAPFeatsMatrix(ValidBoutSyllableIndices{i}(Indices),j))), '*', 'FontSize', 14);
                end
            end
        end
        
        ylabel(BirdParameters(i).SAPFeat_FieldNames{j});
        if (j > 6)
            set(gca, 'XTick', 0:1:length(PlotXLabelString)-1, 'XTickLabel', PlotXLabelString, 'XTickLabelRotation', 45);
        else
            set(gca, 'XTick', 1:1:length(PlotXLabelString), 'XTickLabel', []);
        end        
        if (j < 3)
            title(BirdParameters(i).BirdName);
        end
        axis tight;
        Temp = axis;
        Temp = [-0.5 3.5 0.98*Temp(3) 1.02*Temp(4)];
        axis(Temp);
    end
    set(gcf, 'PaperPositionMode', 'auto');
    print([BirdParameters(i).BirdName, '.SongBoutSyllFeaturesForSyllCategories.png'], '-dpng', '-r300');
end

% Now to plot group data across all birds
for i = 1:length(BirdParameters),
    for j = 1:length(BirdParameters(i).SAPFeat_FieldNames)-1,
        for k = min(ValidBoutSyllableCategories{i}):max(ValidBoutSyllableCategories{i}),
            % Get all the indices for syllables that belong a particular
            % category and part of song bouts
            
            Indices = ValidBoutSyllableIndices{i}(find((ValidBoutSyllableCategories{i}(:) == k) & (ValidBoutSyllableBoutCategories{i} == 1)));
            if (~isempty(Indices))
                MeanFeatureValueAcrossBirds{j}(i,k+1) = mean(BirdParameters(i).SAPFeatsMatrix(Indices,j));
            else
                MeanFeatureValueAcrossBirds{j}(i,k+1) = NaN;
            end
        end
    end
end

disp('Finished');
