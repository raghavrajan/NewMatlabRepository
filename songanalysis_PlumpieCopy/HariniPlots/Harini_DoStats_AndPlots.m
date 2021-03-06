function [OutputStats] = Harini_DoStats_AndPlots(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNot)

% Now to do the stats
% 1. Compare L0 directed and undirected
% 2. For the ones where L0 directed and undirected are significantly
% different, have to do a linear correlation and see if there is a
% correlation or not
% 3. Plot individual L0 directed and undirected comparisons
% 4. Plot average % change across all birds over distance

Colours = 'rgbcmk';

figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [188 208 1600 800]);
for i = 1:length(DirBout_Stats)
    disp(['Bird #', num2str(i), ': ', BirdNames{i}]);
    GroupData{i} = [];
    GLMGroupData{i} = [];
    
    % GroupData has 3 columns, the first one has the data corresponding to
    % the variable that that I'm looking at like TotalINNumber, etc. 2nd
    % column has the index of the distance at which the recordings were
    % done (0 for L0, 1 for L1 ... 6 for UN) and 3rd column has 1 for Dir,
    % 2 for Dir/Undir - DUN - and 3 for Undir
    
    for j = 1:length(DirBout_Stats{i}),
        if (~isempty(DirBout_Stats{i}{j}))
            if (eval(['length(~isnan(DirBout_Stats{i}{j}.', StructString, '(:)))']) >= MinTrialNo)
                GroupData{i} = [GroupData{i}; [eval(['DirBout_Stats{i}{j}.', StructString, '(:)']) ones(size(eval(['DirBout_Stats{i}{j}.', StructString, '(:)'])))*j ones(size(eval(['DirBout_Stats{i}{j}.', StructString, '(:)'])))*1]];
                GLMGroupData{i} = [GLMGroupData{i}; [eval(['DirBout_Stats{i}{j}.', StructString, '(:)']) ones(size(eval(['DirBout_Stats{i}{j}.', StructString, '(:)'])))*j ones(size(eval(['DirBout_Stats{i}{j}.', StructString, '(:)'])))*1]];
            else
                GLMGroupData{i} = [GLMGroupData{i}; [ones(MinTrialNo,1)*NaN ones(MinTrialNo,1)*j ones(MinTrialNo,1)*1]];
            end
        else
            GLMGroupData{i} = [GLMGroupData{i}; [ones(MinTrialNo,1)*NaN ones(MinTrialNo,1)*j ones(MinTrialNo,1)*1]];
        end
    end
    
    for j = 1:length(DirUndirBout_Stats{i}),
        if (~isempty(DirUndirBout_Stats{i}{j}))
            if (eval(['length(~isnan(DirUndirBout_Stats{i}{j}.', StructString, '(:)))']) >= MinTrialNo)
                GroupData{i} = [GroupData{i}; [eval(['DirUndirBout_Stats{i}{j}.', StructString, '(:)']) ones(size(eval(['DirUndirBout_Stats{i}{j}.', StructString, '(:)'])))*j ones(size(eval(['DirUndirBout_Stats{i}{j}.', StructString, '(:)'])))*2]];
                GLMGroupData{i} = [GLMGroupData{i}; [eval(['DirUndirBout_Stats{i}{j}.', StructString, '(:)']) ones(size(eval(['DirUndirBout_Stats{i}{j}.', StructString, '(:)'])))*j ones(size(eval(['DirUndirBout_Stats{i}{j}.', StructString, '(:)'])))*2]];
            else
                GLMGroupData{i} = [GLMGroupData{i}; [ones(MinTrialNo,1)*NaN ones(MinTrialNo,1)*j ones(MinTrialNo,1)*2]];
            end
        else
            GLMGroupData{i} = [GLMGroupData{i}; [ones(MinTrialNo,1)*NaN ones(MinTrialNo,1)*j ones(MinTrialNo,1)*2]];
        end
    end
    
    for j = 1:length(UndirBout_Stats{i}),
        if (~isempty(UndirBout_Stats{i}{j}))
            if (eval(['length(~isnan(UndirBout_Stats{i}{j}.', StructString, '(:)))']) >= MinTrialNo)
                GroupData{i} = [GroupData{i}; [eval(['UndirBout_Stats{i}{j}.', StructString, '(:)']) ones(size(eval(['UndirBout_Stats{i}{j}.', StructString, '(:)'])))*j ones(size(eval(['UndirBout_Stats{i}{j}.', StructString, '(:)'])))*3]];
                GLMGroupData{i} = [GLMGroupData{i}; [eval(['UndirBout_Stats{i}{j}.', StructString, '(:)']) ones(size(eval(['UndirBout_Stats{i}{j}.', StructString, '(:)'])))*j ones(size(eval(['UndirBout_Stats{i}{j}.', StructString, '(:)'])))*3]];
            else
                GLMGroupData{i} = [GLMGroupData{i}; [ones(MinTrialNo,1)*NaN ones(MinTrialNo,1)*j ones(MinTrialNo,1)*3]];
            end
        else
            GLMGroupData{i} = [GLMGroupData{i}; [ones(MinTrialNo,1)*NaN ones(MinTrialNo,1)*j ones(MinTrialNo,1)*3]];
        end
    end

    % Remove NaN values
    NaNIndices = find(isnan(GroupData{i}(:,1)));
    if (~isempty(NaNIndices))
        GroupData{i}(NaNIndices,:) = [];
    end
    
    % Run correlations - pearsons on dir data and undir data
    AllDirSongIndices = find(GroupData{i}(:,3) == 1);
    [r, p] = corrcoef(GroupData{i}(AllDirSongIndices, 1), Distances(GroupData{i}(AllDirSongIndices, 2)));
%    [Coeffs, Stats] = robustfit(GroupData{i}(AllDirSongIndices, 2), GroupData{i}(AllDirSongIndices, 1));
    
    Distance_DirSong_Corr(i,:) = [r(1,2) p(1,2)];
    if (ParametricOrNot == 0)
        [KWP_DirSong{i}, anovatab, KWStats_DirSong{i}] = kruskalwallis(GroupData{i}(AllDirSongIndices, 1), GroupData{i}(AllDirSongIndices, 2), 'off');
    else
        [KWP_DirSong{i}, anovatab, KWStats_DirSong{i}] = anova1(GroupData{i}(AllDirSongIndices, 1), GroupData{i}(AllDirSongIndices, 2), 'off');
    end        
    
    AllUnDirSongIndices = find((GroupData{i}(:,3) == 3) & (GroupData{i}(:,2) < 6));
    [r, p] = corrcoef(GroupData{i}(AllUnDirSongIndices, 1), Distances(GroupData{i}(AllUnDirSongIndices, 2)));
    Distance_UnDirSong_Corr(i,:) = [r(1,2) p(1,2)];
    if (ParametricOrNot == 0)
        [KWP_UnDirSong{i}, anovatab, KWStats_UnDirSong{i}] = kruskalwallis(GroupData{i}(AllUnDirSongIndices, 1), GroupData{i}(AllUnDirSongIndices, 2), 'off');
    else
        [KWP_UnDirSong{i}, anovatab, KWStats_UnDirSong{i}] = anova1(GroupData{i}(AllUnDirSongIndices, 1), GroupData{i}(AllUnDirSongIndices, 2), 'off');
    end
    
    L0Dir_Indices = find((GroupData{i}(:,2) == 1) & (GroupData{i}(:,3) == 1));
    UnDir_Indices = find((GroupData{i}(:,2) == 6) & (GroupData{i}(:,3) == 3));

    if ((length(L0Dir_Indices) >= MinTrialNo) && (length(UnDir_Indices) >= MinTrialNo))
        switch ParametricOrNot
            case 0
                L0Dir_Undir_PValue(i) = kruskalwallis([GroupData{i}(L0Dir_Indices,1); GroupData{i}(UnDir_Indices,1)], [ones(size(L0Dir_Indices(:)))*1; ones(size(UnDir_Indices(:)))*2], 'off');
                % L0Dir_Undir_PValue(i) = signrank(GroupData{i}(L0Dir_Indices,1), GroupData{i}(UnDir_Indices,1));
            case 1
                [HypothesisValue, L0Dir_Undir_PValue(i)] = ttest2(GroupData{i}(L0Dir_Indices,1), GroupData{i}(UnDir_Indices,1));
        end
    else
        L0Dir_Undir_PValue(i) = NaN;
    end
    
    DisplayString = ['P value for ', StructString, ': '];
    if (ParametricOrNot == 0)
        DisplayString = [DisplayString, 'Kruskalwallis test '];
    else
        DisplayString = [DisplayString, 't-test '];
    end
    
    DisplayString = [DisplayString, 'comparing L0 Directed songs and Undirected songs = ', num2str(L0Dir_Undir_PValue(i))];
    
    if (mean(GroupData{i}(L0Dir_Indices,1)) > mean(GroupData{i}(UnDir_Indices,1)))
        DisplayString = [DisplayString, '; Dir > Undir'];
    else
        if (mean(GroupData{i}(L0Dir_Indices,1)) < mean(GroupData{i}(UnDir_Indices,1)))
            DisplayString = [DisplayString, '; Dir < Undir'];
        else
            DisplayString = [DisplayString, '; Dir = Undir'];
        end
    end
    
    disp(DisplayString);
    
    % Now plot the distributions and means for all birds across all
    % distances for directed songs
    subplot(3, ceil(length(DirBout_Stats)/3),i);
    hold on;
    if (~isempty(GroupData{i}))
        Edges = linspace(min(GroupData{i}(:,1)), max(GroupData{i}(:,1)), 5);

        for j = min(GroupData{i}(:,2)):max(GroupData{i}(:,2)),
            DistIndices = find((GroupData{i}(:,2) == j) & (GroupData{i}(:,3) == 1));
            plot(Edges, 100*histc(GroupData{i}(DistIndices,1), Edges)/length(DistIndices), [Colours(j),'o-']);
        end

        axis tight;
        Temp = axis;

        for j = min(GroupData{i}(:,2)):max(GroupData{i}(:,2)),
            DistIndices = find((GroupData{i}(:,2) == j) & (GroupData{i}(:,3) == 1));
            plot(mean(GroupData{i}(DistIndices,1)), (1+j*0.05)*1.05*Temp(4), [Colours(j), 's']);
            plot(mean(GroupData{i}(DistIndices,1)) + [std(GroupData{i}(DistIndices,1))/sqrt(length(DistIndices)) -std(GroupData{i}(DistIndices,1))/sqrt(length(DistIndices))], ones(1,2)*(1+j*0.05)*1.05*Temp(4), Colours(j));
        end

        axis tight;
        xlabel(LabelString);
        ylabel('% of bouts');
        title(BirdNames{i});
    end
    
    if ((~isempty(L0Dir_Indices)) && (~isempty(UnDir_Indices)))
        L0Dir_MeanValues(i) = mean(GroupData{i}(L0Dir_Indices,1));
        L0Dir_MedianValues(i) = median(GroupData{i}(L0Dir_Indices,1));
        L0Dir_MaxValues(i) = max(GroupData{i}(L0Dir_Indices,1));
        UnDir_MeanValues(i) = mean(GroupData{i}(UnDir_Indices,1));
        UnDir_MedianValues(i) = median(GroupData{i}(UnDir_Indices,1));
        UnDir_MaxValues(i) = max(GroupData{i}(UnDir_Indices,1));
    else
        L0Dir_MeanValues(i) = NaN;
        L0Dir_MedianValues(i) = NaN;
        L0Dir_MaxValues(i) = NaN;
        UnDir_MeanValues(i) = NaN;
        UnDir_MedianValues(i) = NaN;
        UnDir_MaxValues(i) = NaN;
    end
    
    TempDataIndices = find(GLMGroupData{i}(:,3) == 2);
    Temp2DataIndices = find(GLMGroupData{i}(:,2) == 6);
    GLMModelDataIndices = setdiff(1:1:size(GroupData{i},1), union(TempDataIndices, Temp2DataIndices));
    % GLMModelDataIndices = setdiff(1:1:size(GLMGroupData{i},1), TempDataIndices);
    switch StructString
        case 'TotalINNumber_500ms'
            DistrName = 'poisson';
        case 'CompleteMotifNumber'
            DistrName = 'poisson';
        case 'FirstMotifDuration'
            DistrName = 'normal';
    end
    
    DistrName = 'normal';
    ContextValues = {'D' 'DUN' 'UN'};
    ContextNumericalValues = [0 2 1];
    
    NumINs = GLMGroupData{i}(GLMModelDataIndices,1);
    Contexts = ContextValues(GLMGroupData{i}(GLMModelDataIndices,3));
    DistanceValues = Distances(GLMGroupData{i}(GLMModelDataIndices,2));
    
    Context_NumericalValues = ContextNumericalValues(GLMGroupData{i}(GLMModelDataIndices, 3));
    
    NumINs = NumINs(:);
    Contexts = Contexts(:);
    DistanceValues = DistanceValues(:);
    
    GLM_DataTable = table(DistanceValues, Contexts, NumINs);
    % mdl = fitglm([Distances(GroupData{i}(GLMModelDataIndices,2))' ContextValues(GroupData{i}(GLMModelDataIndices,3))'], GroupData{i}(GLMModelDataIndices,1), 'interactions', 'CategoricalVars', [2], 'VarNames', {'Distance', 'Context_Values', 'NumINs'}, 'Distribution', DistrName);
    %mdl = fitglm(GLM_DataTable, 'interactions', 'CategoricalVars', [2], 'Distribution', DistrName);
    mdl = fitlm(GLM_DataTable, 'interactions', 'CategoricalVars', [2]);
    
    Temp_Estimates_PValues = mdl.Coefficients{:, {'Estimate', 'pValue'}};
    
%     % if interaction term is not significant, then re-run the model without
%     % interactions
%     
%     if (size(Temp_Estimates_PValues,1) == 4)
%         if (Temp_Estimates_PValues(4,2) >= 0.05)
%             mdl = fitglm(GLM_DataTable, 'linear', 'CategoricalVars', [2], 'Distribution', DistrName);
%             mdl2 = fitglm([DistanceValues(:) Context_NumericalValues(:)], NumINs, 'linear', 'CategoricalVars', [2], 'Distribution', DistrName);
%             mdl3 = fitglm(zscore([DistanceValues(:) Context_NumericalValues(:)]), (NumINs(:) - nanmean(NumINs))/nanstd(NumINs), 'linear', 'CategoricalVars', [2], 'Distribution', DistrName);
%         else
%             mdl2 = fitglm([DistanceValues(:) Context_NumericalValues(:)], NumINs, 'interactions', 'CategoricalVars', [2], 'Distribution', DistrName);
%             mdl3 = fitglm(zscore([DistanceValues(:) Context_NumericalValues(:)]), (NumINs(:) - nanmean(NumINs))/nanstd(NumINs), 'interactions', 'CategoricalVars', [2], 'Distribution', DistrName);
%         end    
%     else
%         mdl2 = fitglm([DistanceValues(:) Context_NumericalValues(:)], NumINs, 'interactions', 'CategoricalVars', [2], 'Distribution', DistrName);
%         mdl3 = fitglm(zscore([DistanceValues(:) Context_NumericalValues(:)]), (NumINs(:) - nanmean(NumINs))/nanstd(NumINs), 'interactions', 'CategoricalVars', [2], 'Distribution', DistrName);
%     end
    
    %mdl3 = fitglm(zscore([DistanceValues(:) Context_NumericalValues(:)]), (NumINs(:) - nanmean(NumINs))/nanstd(NumINs), 'interactions', 'CategoricalVars', [2], 'Distribution', DistrName);
    mdl3 = fitlm(zscore([DistanceValues(:) Context_NumericalValues(:)]), (NumINs(:) - nanmean(NumINs))/nanstd(NumINs), 'interactions', 'CategoricalVars', [2]);
    Temp_Estimates_PValues = mdl3.Coefficients{:, {'Estimate', 'pValue'}};
    
    if (size(Temp_Estimates_PValues,1) == 4)
        GLMFit_Intercept(i,:) = Temp_Estimates_PValues(1,:);
        GLMFit_Distance(i,:) = Temp_Estimates_PValues(2,:);
        GLMFit_Context(i,:) = Temp_Estimates_PValues(3,:);
        GLMFit_DistanceContext_Interactions(i,:) = Temp_Estimates_PValues(4,:);
    else
        if (size(Temp_Estimates_PValues,1) == 3)
            GLMFit_Intercept(i,:) = Temp_Estimates_PValues(1,:);
            GLMFit_Distance(i,:) = Temp_Estimates_PValues(2,:);
            GLMFit_Context(i,:) = Temp_Estimates_PValues(3,:);
            GLMFit_DistanceContext_Interactions(i,:) = ones(1,2)*NaN;
        else
            GLMFit_Intercept(i,:) = ones(1,2)*NaN;
            GLMFit_Distance(i,:) = ones(1,2)*NaN;
            GLMFit_Context(i,:) = ones(1,2)*NaN;
            GLMFit_DistanceContext_Interactions(i,:) = ones(1,2)*NaN;
        end
    end
    disp(mdl);
    disp(mdl3);
end

disp('Distance directed song correlations: ');
disp(Distance_DirSong_Corr);

disp('Distance undirected song correlations: ');
disp(Distance_UnDirSong_Corr);

disp('Group data comparisons');
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [131 558 1100 400]);
p = panel();
p.pack('h', {1/3 1/3 1/3});
p.de.margin = 10;

p(1).select();
hold on;
L0DirBar = bar(1, nanmean(L0Dir_MeanValues));
set(L0DirBar, 'FaceColor', 'none', 'LineWidth', 1);
errorbar(1, nanmean(L0Dir_MeanValues), nanstd(L0Dir_MeanValues)/sqrt(length(find(~isnan(L0Dir_MeanValues)))), 'k', 'LineWidth', 1);

UnDirBar = bar(2, nanmean(UnDir_MeanValues));
set(UnDirBar, 'FaceColor', 'none', 'LineWidth', 1);
errorbar(2, nanmean(UnDir_MeanValues), nanstd(UnDir_MeanValues)/sqrt(length(find(~isnan(UnDir_MeanValues)))), 'k', 'LineWidth', 1);

for i = 1:length(UnDir_MeanValues),
    if (L0Dir_Undir_PValue(i) < 0.05)
        plot([1.1 1.9], [L0Dir_MeanValues(i) UnDir_MeanValues(i)], 'ko-', 'MarkerFaceColor', 'k');
    else
        if (~isnan(L0Dir_Undir_PValue(i)))
            plot([1.1 1.9], [L0Dir_MeanValues(i) UnDir_MeanValues(i)], 'ko-');
        end
    end
end
set(gca, 'XTick', [1 2], 'XTickLabel', {'L0 Dir' 'Undir'});
axis auto;
ylabel(['Mean of ', LabelString]);

% Stats
Group_L0Dir_UnDir_PValue.Mean = signrank(L0Dir_MeanValues(find(~isnan(L0Dir_MeanValues))), UnDir_MeanValues(find(~isnan(L0Dir_MeanValues))));
sigstar({[1 2]}, Group_L0Dir_UnDir_PValue.Mean);
if (nanmean(L0Dir_MeanValues) > nanmean(UnDir_MeanValues))
    disp(['Sign rank test comparing group means for L0 directed and Undirected, p = ', num2str(Group_L0Dir_UnDir_PValue.Mean), '; Dir > Undir']);
else
    if (nanmean(L0Dir_MeanValues) < nanmean(UnDir_MeanValues))
        disp(['Sign rank test comparing group means for L0 directed and Undirected, p = ', num2str(Group_L0Dir_UnDir_PValue.Mean), '; Dir < Undir']);
    else
        disp(['Sign rank test comparing group means for L0 directed and Undirected, p = ', num2str(Group_L0Dir_UnDir_PValue.Mean), '; Dir = Undir']);
    end
end

p(2).select();
hold on;
L0DirBar = bar(1, nanmean(L0Dir_MedianValues));
set(L0DirBar, 'FaceColor', 'none', 'LineWidth', 1);
errorbar(1, nanmean(L0Dir_MedianValues), nanstd(L0Dir_MedianValues)/sqrt(length(find(~isnan(L0Dir_MedianValues)))), 'k', 'LineWidth', 1);

UnDirBar = bar(2, nanmean(UnDir_MedianValues));
set(UnDirBar, 'FaceColor', 'none', 'LineWidth', 1);
errorbar(2, nanmean(UnDir_MedianValues), nanstd(UnDir_MedianValues)/sqrt(length(find(~isnan(UnDir_MedianValues)))), 'k', 'LineWidth', 1);

for i = 1:length(UnDir_MedianValues),
    if (L0Dir_Undir_PValue(i) < 0.05)
        plot([1.1 1.9], [L0Dir_MedianValues(i) UnDir_MedianValues(i)] + ones(1,2)*rand/5, 'ko-', 'MarkerFaceColor', 'k');
    else
        if (~isnan(L0Dir_Undir_PValue(i)))
            plot([1.1 1.9], [L0Dir_MedianValues(i) UnDir_MedianValues(i)] + ones(1,2)*rand/5, 'ko-');
        end
    end
end
set(gca, 'XTick', [1 2], 'XTickLabel', {'L0 Dir' 'Undir'});
axis auto;
ylabel(['Median of ', LabelString]);
% Stats
Group_L0Dir_UnDir_PValue.Median = signrank(L0Dir_MedianValues(find(~isnan(L0Dir_MedianValues))), UnDir_MedianValues(find(~isnan(L0Dir_MedianValues))));
sigstar({[1 2]}, Group_L0Dir_UnDir_PValue.Median);
if (nanmean(L0Dir_MedianValues) > nanmean(UnDir_MeanValues))
    disp(['Sign rank test comparing group medians for L0 directed and Undirected, p = ', num2str(Group_L0Dir_UnDir_PValue.Median), '; Dir > Undir']);
else
    if (nanmean(L0Dir_MedianValues) < nanmean(UnDir_MedianValues))
        disp(['Sign rank test comparing group medians for L0 directed and Undirected, p = ', num2str(Group_L0Dir_UnDir_PValue.Median), '; Dir < Undir']);
    else
        disp(['Sign rank test comparing group medians for L0 directed and Undirected, p = ', num2str(Group_L0Dir_UnDir_PValue.Median), '; Dir = Undir']);
    end
end

p(3).select();
hold on;
L0DirBar = bar(1, nanmean(L0Dir_MaxValues));
set(L0DirBar, 'FaceColor', 'none', 'LineWidth', 1);
errorbar(1, nanmean(L0Dir_MaxValues), nanstd(L0Dir_MaxValues)/sqrt(length(find(~isnan(L0Dir_MaxValues)))), 'k', 'LineWidth', 1);

UnDirBar = bar(2, nanmean(UnDir_MaxValues));
set(UnDirBar, 'FaceColor', 'none', 'LineWidth', 1);
errorbar(2, nanmean(UnDir_MaxValues), nanstd(UnDir_MaxValues)/sqrt(length(find(~isnan(UnDir_MaxValues)))), 'k', 'LineWidth', 1);

for i = 1:length(UnDir_MaxValues),
    if (L0Dir_Undir_PValue(i) < 0.05)
        plot([1.1 1.9], [L0Dir_MaxValues(i) UnDir_MaxValues(i)] + ones(1,2)*rand/5, 'ko-', 'MarkerFaceColor', 'k');
    else
        if (~isnan(L0Dir_Undir_PValue(i)))
            plot([1.1 1.9], [L0Dir_MaxValues(i) UnDir_MaxValues(i)] + ones(1,2)*rand/5, 'ko-');
        end
    end
end
set(gca, 'XTick', [1 2], 'XTickLabel', {'L0 Dir' 'Undir'});
axis auto;
ylabel(['Max of ', LabelString]);
% Stats
Group_L0Dir_UnDir_PValue.Max = signrank(L0Dir_MaxValues(find(~isnan(L0Dir_MaxValues))), UnDir_MaxValues(find(~isnan(L0Dir_MaxValues))));
sigstar({[1 2]}, Group_L0Dir_UnDir_PValue.Max);
if (nanmean(L0Dir_MaxValues) > nanmean(UnDir_MaxValues))
    disp(['Sign rank test comparing group max values for L0 directed and Undirected, p = ', num2str(Group_L0Dir_UnDir_PValue.Max), '; Dir > Undir']);
else
    if (nanmean(L0Dir_MaxValues) < nanmean(UnDir_MaxValues))
        disp(['Sign rank test comparing group max values for L0 directed and Undirected, p = ', num2str(Group_L0Dir_UnDir_PValue.Max), '; Dir < Undir']);
    else
        disp(['Sign rank test comparing group max values for L0 directed and Undirected, p = ', num2str(Group_L0Dir_UnDir_PValue.Max), '; Dir = Undir']);
    end
end


% Now to plot the mean and std across distances for D, D/UN and UN for all
% birds together as percent change relative to undirected.
% This is done only for birds where the means are already different between
% L0 directed and undirected as measured by the earlier p-value

GroupDataFig = figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [800 100 500 800]);
p = panel();
p.pack({1/3 1/3 1/3});
p.de.margin = 15;
p(1).select();
hold on;
p(2).select();
hold on;
p(3).select();
hold on;

BirdsWithSigDiffs = find(L0Dir_Undir_PValue < 0.05);
BirdsWithoutSigDiffs = (find((L0Dir_Undir_PValue >= 0.05) | (isnan(L0Dir_Undir_PValue))));

DistanceDirMeans = [];
DistanceDirUnDirMeans = [];
DistanceUnDirMeans = [];
Symbols = 'so^+hp><*dxv.so';

Index = 0;

DirUnDirComparisonFig = figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [188 208 1600 800]);
for i = 1:length(BirdNames),
    Index = i;
    AllDirSongIndices = find(GroupData{i}(:,3) == 1);
    
    % First for directed songs
    for j = 1:6,
        DirSongIndices = find((GroupData{i}(:,3) == 1) & (GroupData{i}(:,2) == j));
        DirSong_NumTrials(i,j) = length(DirSongIndices);
        if (length(DirSongIndices) >= MinTrialNo)
            DistanceDirMeans(Index, j) = mean(GroupData{i}(DirSongIndices,1));
            DistanceDirSEMs(Index, j) = std(GroupData{i}(DirSongIndices,1))/sqrt(length(DirSongIndices));
            DistanceDirCIs{Index}(j,:) = bootci(10000,@mean, GroupData{i}(DirSongIndices,1));
        else
            DistanceDirMeans(Index, j) = NaN;
            DistanceDirSEMs(Index,j) = NaN;
            DistanceDirCIs{Index}(j,:) = ones(1,2)*NaN;
        end
        
        DirUnDirSongIndices = find((GroupData{i}(:,3) == 2) & (GroupData{i}(:,2) == j));
        DirUnDirSong_NumTrials(i,j) = length(DirUnDirSongIndices);
        if (length(DirUnDirSongIndices) >= MinTrialNo)
            DistanceDirUnDirMeans(Index, j) = mean(GroupData{i}(DirUnDirSongIndices,1));
            DistanceDirUnDirSEMs(Index, j) = std(GroupData{i}(DirUnDirSongIndices,1))/sqrt(length(DirUnDirSongIndices));
            DistanceDirUnDirCIs{Index}(j,:) = bootci(10000, @mean, GroupData{i}(DirUnDirSongIndices,1));
        else
            DistanceDirUnDirMeans(Index, j) = NaN;
            DistanceDirUnDirSEMs(Index, j) = NaN;
            DistanceDirUnDirCIs{Index}(j,:) = ones(1,2)*NaN;
        end
        
        UnDirSongIndices = find((GroupData{i}(:,3) == 3) & (GroupData{i}(:,2) == j));
        UnDirSong_NumTrials(i,j) = length(UnDirSongIndices);
        if (length(UnDirSongIndices) >= MinTrialNo)
            DistanceUnDirMeans(Index, j) = mean(GroupData{i}(UnDirSongIndices,1));
            DistanceUnDirSEMs(Index, j) = std(GroupData{i}(UnDirSongIndices,1))/sqrt(length(UnDirSongIndices));
            DistanceUnDirCIs{Index}(j,:) = bootci(10000, @mean, GroupData{i}(UnDirSongIndices,1));
        else
            DistanceUnDirMeans(Index, j) = NaN;
            DistanceUnDirSEMs(Index, j) = NaN;
            DistanceUnDirCIs{Index}(j,:) = ones(1,2)*NaN;
        end
    end
    
    figure(DirUnDirComparisonFig);
    subplot(3, ceil(length(DirBout_Stats)/3),Index);
    hold on;
    plot(Distances, DistanceDirMeans(Index,:), 'rs-', 'LineWidth', 1);
    plot(repmat(Distances(:), 1, 2)', DistanceDirCIs{Index}', 'r')
    
%     plot(Distances, DistanceDirUnDirMeans(Index,:), 'ks-', 'LineWidth', 1);
%     plot(repmat(Distances(:), 1, 2)', DistanceDirUnDirCIs{Index}', 'k')
%     
    plot(Distances, DistanceUnDirMeans(Index,:), 'bs-', 'LineWidth', 1);
    plot(repmat(Distances(:), 1, 2)', DistanceUnDirCIs{Index}', 'b');
    plot(repmat(Distances(:),1,2), repmat(DistanceUnDirCIs{Index}(end,:), length(Distances), 1), 'k--');
    
    xlabel('Distance from female (cm)');
    axis tight;
    Temp = axis;
    Temp = [-5 250 0.95*Temp(3) 1.05*Temp(4)];
    axis(Temp);
    ylabel(LabelString);
    title({[BirdNames{i}, ': r=', num2str(Distance_DirSong_Corr(i,1)), '; p=', num2str(Distance_DirSong_Corr(i,2)), ';p=', num2str(KWP_DirSong{i})]; [BirdNames{i}, ': r=', num2str(Distance_UnDirSong_Corr(i,1)), '; p=', num2str(Distance_UnDirSong_Corr(i,2)), ';p=', num2str(KWP_UnDirSong{i})]}, 'Color', 'k', 'FontSize', 8, 'FontWeight', 'normal');
    % legend('D', 'DUN', 'UN');
    
%     % Set the value for directed songs and dir/undir songs for the
%     % undirected to be equal to the undirected song means
%        
%     DistanceDirMeans(Index,end) = DistanceUnDirMeans(Index,end);
%     DistanceDirUnDirMeans(Index,end) = DistanceUnDirMeans(Index,end);
%     
%     
%     DistanceDirMeans(Index,:) = 100*DistanceDirMeans(Index,:)/DistanceDirMeans(Index,end);
%     DistanceDirUnDirMeans(Index,:) = 100*DistanceDirUnDirMeans(Index,:)/DistanceDirUnDirMeans(Index,end);
%     DistanceUnDirMeans(Index,:) = 100*DistanceUnDirMeans(Index,:)/DistanceUnDirMeans(Index,end);
    
    if (~isempty(find(BirdsWithSigDiffs == i)))
        figure(GroupDataFig);
        if (L0Dir_MeanValues(i) > UnDir_MeanValues(i))
            p(1).select();
            plot(Distances, DistanceDirMeans(Index,:), ['r', Symbols(i), '-']);

            p(2).select();
            plot(Distances, DistanceDirUnDirMeans(Index,:), ['r', Symbols(i), '-']);

            p(3).select();
            plot(Distances, DistanceUnDirMeans(Index,:), ['r', Symbols(i), '-']);
        else
            if (L0Dir_MeanValues(i) < UnDir_MeanValues(i))
                p(1).select();
                plot(Distances, DistanceDirMeans(Index,:), ['b', Symbols(i), '-']);

                p(2).select();
                plot(Distances, DistanceDirUnDirMeans(Index,:), ['b', Symbols(i), '-']);

                p(3).select();
                plot(Distances, DistanceUnDirMeans(Index,:), ['b', Symbols(i), '-']);
            else
                p(1).select();
                title('Directed songs');
                plot(Distances, DistanceDirMeans(Index,:), ['k', Symbols(i), '-']);

                p(2).select();
                title('DUN songs');
                plot(Distances, DistanceDirUnDirMeans(Index,:), ['k', Symbols(i), '-']);

                p(3).select();
                title('Undir songs');
                plot(Distances, DistanceUnDirMeans(Index,:), ['k', Symbols(i), '-']);
            end
        end
    end
end

% for i = (BirdsWithoutSigDiffs(:))',
%     Index = Index + 1;
%     % First for directed songs
%     for j = 1:6,
%         DirSongIndices = find((GroupData{i}(:,3) == 1) & (GroupData{i}(:,2) == j));
%         DirSong_NumTrials(i,j) = length(DirSongIndices);
%         if (length(DirSongIndices) >= MinTrialNo)
%             DistanceDirMeans(Index, j) = mean(GroupData{i}(DirSongIndices,1));
%             DistanceDirSEMs(Index, j) = std(GroupData{i}(DirSongIndices,1))/sqrt(length(DirSongIndices));
%             DistanceDirCIs{Index}(j,:) = bootci(10000,@mean, GroupData{i}(DirSongIndices,1));
%         else
%             DistanceDirMeans(Index, j) = NaN;
%             DistanceDirSEMs(Index,j) = NaN;
%             DistanceDirCIs{Index}(j,:) = ones(1,2)*NaN;
%         end
% 
%         DirUnDirSongIndices = find((GroupData{i}(:,3) == 2) & (GroupData{i}(:,2) == j));
%         DirUnDirSong_NumTrials(i,j) = length(DirUnDirSongIndices);
%         if (length(DirUnDirSongIndices) >= MinTrialNo)
%             DistanceDirUnDirMeans(Index, j) = mean(GroupData{i}(DirUnDirSongIndices,1));
%             DistanceDirUnDirSEMs(Index, j) = std(GroupData{i}(DirUnDirSongIndices,1))/sqrt(length(DirUnDirSongIndices));
%             DistanceDirUnDirCIs{Index}(j,:) = bootci(10000, @mean, GroupData{i}(DirUnDirSongIndices,1));
%         else
%             DistanceDirUnDirMeans(Index, j) = NaN;
%             DistanceDirUnDirSEMs(Index, j) = NaN;
%             DistanceDirUnDirCIs{Index}(j,:) = ones(1,2)*NaN;
%         end
% 
%         UnDirSongIndices = find((GroupData{i}(:,3) == 3) & (GroupData{i}(:,2) == j));
%         UnDirSong_NumTrials(i,j) = length(UnDirSongIndices);
%         if (length(UnDirSongIndices) >= MinTrialNo)
%             DistanceUnDirMeans(Index, j) = mean(GroupData{i}(UnDirSongIndices,1));
%             DistanceUnDirSEMs(Index, j) = std(GroupData{i}(UnDirSongIndices,1))/sqrt(length(UnDirSongIndices));
%             DistanceUnDirCIs{Index}(j,:) = bootci(10000, @mean, GroupData{i}(UnDirSongIndices,1));
%         else
%             DistanceUnDirMeans(Index, j) = NaN;
%             DistanceUnDirSEMs(Index, j) = NaN;
%             DistanceUnDirCIs{Index}(j,:) = ones(1,2)*NaN;
%         end
%     end
% 
%     % Set the value for directed songs and dir/undir songs for the
%     % undirected to be equal to the undirected song means
%        
%     DistanceDirMeans(Index,end) = DistanceUnDirMeans(Index,end);
%     DistanceDirUnDirMeans(Index,end) = DistanceUnDirMeans(Index,end);
%     
%     
%     DistanceDirMeans(Index,:) = 100*DistanceDirMeans(Index,:)/DistanceDirMeans(Index,end);
%     DistanceDirUnDirMeans(Index,:) = 100*DistanceDirUnDirMeans(Index,:)/DistanceDirUnDirMeans(Index,end);
%     DistanceUnDirMeans(Index,:) = 100*DistanceUnDirMeans(Index,:)/DistanceUnDirMeans(Index,end);
%     
%     figure(DirUnDirComparisonFig);
%     subplot(3, ceil(length(DirBout_Stats)/3),Index);
%     hold on;
%     plot(Distances, DistanceDirMeans(Index,:), 'rs-', 'LineWidth', 1);
%     plot(repmat(Distances(:), 1, 2)', DistanceDirCIs{Index}', 'r')
%     
% %     plot(Distances, DistanceDirUnDirMeans(Index,:), 'ks-', 'LineWidth', 1);
% %     plot(repmat(Distances(:), 1, 2)', DistanceDirUnDirCIs{Index}', 'k')
% %     
%     plot(Distances, DistanceUnDirMeans(Index,:), 'bs-', 'LineWidth', 1);
%     plot(repmat(Distances(:), 1, 2)', DistanceUnDirCIs{Index}', 'b');
%     plot(repmat(Distances(:),1,2), repmat(DistanceUnDirCIs{Index}(end,:), length(Distances), 1), 'k--');
%     
%     xlabel('Distance from female (cm)');
%     axis tight;
%     Temp = axis;
%     Temp = [-0.5 250 0.95*Temp(3) 1.05*Temp(4)];
%     axis(Temp);
%     ylabel(LabelString);
%     title({[BirdNames{i}, ': r=', num2str(Distance_DirSong_Corr(i,1)), '; p=', num2str(Distance_DirSong_Corr(i,2)), ';p=', num2str(KWP_DirSong{i})]; [BirdNames{i}, ': r=', num2str(Distance_UnDirSong_Corr(i,1)), '; p=', num2str(Distance_UnDirSong_Corr(i,2)), ';p=', num2str(KWP_UnDirSong{i})]}, 'Color', 'k', 'FontSize', 8, 'FontWeight', 'normal');
%     % legend('D', 'UN');
% end

figure(GroupDataFig);
for i = 1:3,
    p(i).select();
    if (i == 2)
        ylabel(['% normalized change in ', LabelString]);
    end
    if (i == 3)
        xlabel('Distance to female (cm)');
    end
end

OutputStats.DistanceDirMeans = DistanceDirMeans;
OutputStats.DistanceUnDirMeans = DistanceUnDirMeans;

OutputStats.DistanceDirCIs = DistanceDirCIs;
OutputStats.DistanceUnDirCIs = DistanceUnDirCIs;

OutputStats.Distance_DirSong_Corr = Distance_DirSong_Corr;
OutputStats.Distance_UnDirSong_Corr = Distance_UnDirSong_Corr;

OutputStats.KWP_DirSong = KWP_DirSong;
OutputStats.KWP_UnDirSong = KWP_UnDirSong;

OutputStats.L0Dir_Undir_PValue = L0Dir_Undir_PValue;

OutputStats.GLMFit_Intercept = GLMFit_Intercept;
OutputStats.GLMFit_DistanceContext_Interactions = GLMFit_DistanceContext_Interactions;
OutputStats.GLMFit_Distance = GLMFit_Distance;
OutputStats.GLMFit_Context = GLMFit_Context;

disp('Done with stats and plots');
