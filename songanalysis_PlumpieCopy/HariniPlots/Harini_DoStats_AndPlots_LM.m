function [OutputStats] = Harini_DoStats_AndPlots_LM(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNot)

% Now to do the stats
% 1. Compare L0 directed and undirected
% 2. For the ones where L0 directed and undirected are significantly
% different, have to do a linear correlation and see if there is a
% correlation or not
% 3. Plot individual L0 directed and undirected comparisons
% 4. Plot average % change across all birds over distance

Colours = 'rgbcmk';

for i = 1:length(DirBout_Stats)
    disp(['Bird #', num2str(i), ': ', BirdNames{i}]);
    GroupData{i} = [];
    LMGroupData{i} = [];
    
    % GroupData has 3 columns, the first one has the data corresponding to
    % the variable that that I'm looking at like TotalINNumber, etc. 2nd
    % column has the index of the distance at which the recordings were
    % done (0 for L0, 1 for L1 ... 6 for UN) and 3rd column has 1 for Dir,
    % 2 for Dir/Undir - DUN - and 3 for Undir
    
    for j = 1:length(DirBout_Stats{i}),
        if (~isempty(DirBout_Stats{i}{j}))
            if (eval(['length(~isnan(DirBout_Stats{i}{j}.', StructString, '(:)))']) >= MinTrialNo)
                GroupData{i} = [GroupData{i}; [eval(['DirBout_Stats{i}{j}.', StructString, '(:)']) ones(size(eval(['DirBout_Stats{i}{j}.', StructString, '(:)'])))*j ones(size(eval(['DirBout_Stats{i}{j}.', StructString, '(:)'])))*1]];
                LMGroupData{i} = [LMGroupData{i}; [eval(['DirBout_Stats{i}{j}.', StructString, '(:)']) ones(size(eval(['DirBout_Stats{i}{j}.', StructString, '(:)'])))*j ones(size(eval(['DirBout_Stats{i}{j}.', StructString, '(:)'])))*1]];
            else
                LMGroupData{i} = [LMGroupData{i}; [ones(MinTrialNo,1)*NaN ones(MinTrialNo,1)*j ones(MinTrialNo,1)*1]];
            end
        else
            LMGroupData{i} = [LMGroupData{i}; [ones(MinTrialNo,1)*NaN ones(MinTrialNo,1)*j ones(MinTrialNo,1)*1]];
        end
    end
    
    for j = 1:length(DirUndirBout_Stats{i}),
        if (~isempty(DirUndirBout_Stats{i}{j}))
            if (eval(['length(~isnan(DirUndirBout_Stats{i}{j}.', StructString, '(:)))']) >= MinTrialNo)
                GroupData{i} = [GroupData{i}; [eval(['DirUndirBout_Stats{i}{j}.', StructString, '(:)']) ones(size(eval(['DirUndirBout_Stats{i}{j}.', StructString, '(:)'])))*j ones(size(eval(['DirUndirBout_Stats{i}{j}.', StructString, '(:)'])))*2]];
                LMGroupData{i} = [LMGroupData{i}; [eval(['DirUndirBout_Stats{i}{j}.', StructString, '(:)']) ones(size(eval(['DirUndirBout_Stats{i}{j}.', StructString, '(:)'])))*j ones(size(eval(['DirUndirBout_Stats{i}{j}.', StructString, '(:)'])))*2]];
            else
                LMGroupData{i} = [LMGroupData{i}; [ones(MinTrialNo,1)*NaN ones(MinTrialNo,1)*j ones(MinTrialNo,1)*2]];
            end
        else
            LMGroupData{i} = [LMGroupData{i}; [ones(MinTrialNo,1)*NaN ones(MinTrialNo,1)*j ones(MinTrialNo,1)*2]];
        end
    end
    
    for j = 1:length(UndirBout_Stats{i}),
        if (~isempty(UndirBout_Stats{i}{j}))
            if (eval(['length(~isnan(UndirBout_Stats{i}{j}.', StructString, '(:)))']) >= MinTrialNo)
                GroupData{i} = [GroupData{i}; [eval(['UndirBout_Stats{i}{j}.', StructString, '(:)']) ones(size(eval(['UndirBout_Stats{i}{j}.', StructString, '(:)'])))*j ones(size(eval(['UndirBout_Stats{i}{j}.', StructString, '(:)'])))*3]];
                LMGroupData{i} = [LMGroupData{i}; [eval(['UndirBout_Stats{i}{j}.', StructString, '(:)']) ones(size(eval(['UndirBout_Stats{i}{j}.', StructString, '(:)'])))*j ones(size(eval(['UndirBout_Stats{i}{j}.', StructString, '(:)'])))*3]];
            else
                LMGroupData{i} = [LMGroupData{i}; [ones(MinTrialNo,1)*NaN ones(MinTrialNo,1)*j ones(MinTrialNo,1)*3]];
            end
        else
            LMGroupData{i} = [LMGroupData{i}; [ones(MinTrialNo,1)*NaN ones(MinTrialNo,1)*j ones(MinTrialNo,1)*3]];
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
        [KWP_DirSong{i}, anovatab, KWStats_DirSong{i}] = anova1(GroupData{i}(AllDirSongIndices, 1), GroupData{i}(AllDirSongIndices, 2), 'off');
    else
        [KWP_DirSong{i}, anovatab, KWStats_DirSong{i}] = anova1(GroupData{i}(AllDirSongIndices, 1), GroupData{i}(AllDirSongIndices, 2), 'off');
    end        
    
    AllUnDirSongIndices = find((GroupData{i}(:,3) == 3) & (GroupData{i}(:,2) < 6));
    [r, p] = corrcoef(GroupData{i}(AllUnDirSongIndices, 1), Distances(GroupData{i}(AllUnDirSongIndices, 2)));
    Distance_UnDirSong_Corr(i,:) = [r(1,2) p(1,2)];
    if (ParametricOrNot == 0)
        [KWP_UnDirSong{i}, anovatab, KWStats_UnDirSong{i}] = anova1(GroupData{i}(AllUnDirSongIndices, 1), GroupData{i}(AllUnDirSongIndices, 2), 'off');
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
    
    TempDataIndices = find(LMGroupData{i}(:,3) == 2);
    Temp2DataIndices = find(LMGroupData{i}(:,2) == 6);
    LMModelDataIndices = setdiff(1:1:size(GroupData{i},1), union(TempDataIndices, Temp2DataIndices));
    
    ContextValues = {'D' 'DUN' 'UN'};
    ContextNumericalValues = [0 2 1];
    
    NumINs = LMGroupData{i}(LMModelDataIndices,1);
    Contexts = ContextValues(LMGroupData{i}(LMModelDataIndices,3));
    DistanceValues = Distances(LMGroupData{i}(LMModelDataIndices,2));
    
    Context_NumericalValues = ContextNumericalValues(LMGroupData{i}(LMModelDataIndices, 3));
    
    NumINs = NumINs(:);
    Contexts = Contexts(:);
    DistanceValues = DistanceValues(:);
    
    LM_DataTable = table(DistanceValues, Contexts, NumINs);

    mdl = fitlm(LM_DataTable, 'interactions', 'CategoricalVars', [2]);
    
    Temp_Estimates_PValues = mdl.Coefficients{:, {'Estimate', 'pValue'}};
    
    mdl3 = fitlm(([zscore(DistanceValues(:)) Context_NumericalValues(:)]), (NumINs(:) - nanmean(NumINs))/nanstd(NumINs), 'interactions', 'CategoricalVars', [2]);
    Temp_Estimates_PValues = mdl.Coefficients{:, {'Estimate', 'pValue'}};
    if (size(Temp_Estimates_PValues,1) < 4)
        LMFit_ModelPValue(i) = NaN;
    else
        Temp_Model_Summary = anova(mdl, 'summary');
        LMFit_ModelPValue(i) = table2array(Temp_Model_Summary(2,5));
    end
    
    if (size(Temp_Estimates_PValues,1) == 4)
        LMFit_Intercept(i,:) = Temp_Estimates_PValues(1,:);
        LMFit_Distance(i,:) = Temp_Estimates_PValues(2,:);
        LMFit_Context(i,:) = Temp_Estimates_PValues(3,:);
        LMFit_DistanceContext_Interactions(i,:) = Temp_Estimates_PValues(4,:);
    else
        if (size(Temp_Estimates_PValues,1) == 3)
            LMFit_Intercept(i,:) = Temp_Estimates_PValues(1,:);
            LMFit_Distance(i,:) = Temp_Estimates_PValues(2,:);
            LMFit_Context(i,:) = Temp_Estimates_PValues(3,:);
            LMFit_DistanceContext_Interactions(i,:) = ones(1,2)*NaN;
        else
            LMFit_Intercept(i,:) = ones(1,2)*NaN;
            LMFit_Distance(i,:) = ones(1,2)*NaN;
            LMFit_Context(i,:) = ones(1,2)*NaN;
            LMFit_DistanceContext_Interactions(i,:) = ones(1,2)*NaN;
        end
    end
    
    Temp_Estimates_PValues = mdl3.Coefficients{:, {'Estimate', 'pValue'}};
        if (size(Temp_Estimates_PValues,1) == 4)
        LMFit_ZScore_Intercept(i,:) = Temp_Estimates_PValues(1,:);
        LMFit_ZScore_Distance(i,:) = Temp_Estimates_PValues(2,:);
        LMFit_ZScore_Context(i,:) = Temp_Estimates_PValues(3,:);
        LMFit_ZScore_DistanceContext_Interactions(i,:) = Temp_Estimates_PValues(4,:);
    else
        if (size(Temp_Estimates_PValues,1) == 3)
            LMFit_ZScore_Intercept(i,:) = Temp_Estimates_PValues(1,:);
            LMFit_ZScore_Distance(i,:) = Temp_Estimates_PValues(2,:);
            LMFit_ZScore_Context(i,:) = Temp_Estimates_PValues(3,:);
            LMFit_ZScore_DistanceContext_Interactions(i,:) = ones(1,2)*NaN;
        else
            LMFit_ZScore_Intercept(i,:) = ones(1,2)*NaN;
            LMFit_ZScore_Distance(i,:) = ones(1,2)*NaN;
            LMFit_ZScore_Context(i,:) = ones(1,2)*NaN;
            LMFit_ZScore_DistanceContext_Interactions(i,:) = ones(1,2)*NaN;
        end
    end
    
    Models{i}.mdl = mdl;
    Models{i}.ZScoremdl = mdl3;
    % disp(mdl);
    % disp(mdl3);
end

disp('Distance directed song correlations: ');
disp(Distance_DirSong_Corr);

disp('Distance undirected song correlations: ');
disp(Distance_UnDirSong_Corr);

disp('Group data comparisons');

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

BirdsWithSigDiffs = find(L0Dir_Undir_PValue < 0.05);
BirdsWithoutSigDiffs = (find((L0Dir_Undir_PValue >= 0.05) | (isnan(L0Dir_Undir_PValue))));

DistanceDirMeans = [];
DistanceDirUnDirMeans = [];
DistanceUnDirMeans = [];
Symbols = 'so^+hp><*dxv.so';

Index = 0;

% DirUnDirComparisonFig = figure;
% set(gcf, 'Color', 'w');
% set(gcf, 'Position', [188 208 1600 800]);
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
            % DistanceDirCIs{Index}(j,:) = bootci(10000,@mean, GroupData{i}(DirSongIndices,1));
        else
            DistanceDirMeans(Index, j) = NaN;
            DistanceDirSEMs(Index,j) = NaN;
            % DistanceDirCIs{Index}(j,:) = ones(1,2)*NaN;
        end
        
        DirUnDirSongIndices = find((GroupData{i}(:,3) == 2) & (GroupData{i}(:,2) == j));
        DirUnDirSong_NumTrials(i,j) = length(DirUnDirSongIndices);
        if (length(DirUnDirSongIndices) >= MinTrialNo)
            DistanceDirUnDirMeans(Index, j) = mean(GroupData{i}(DirUnDirSongIndices,1));
            DistanceDirUnDirSEMs(Index, j) = std(GroupData{i}(DirUnDirSongIndices,1))/sqrt(length(DirUnDirSongIndices));
            % DistanceDirUnDirCIs{Index}(j,:) = bootci(10000, @mean, GroupData{i}(DirUnDirSongIndices,1));
        else
            DistanceDirUnDirMeans(Index, j) = NaN;
            DistanceDirUnDirSEMs(Index, j) = NaN;
            % DistanceDirUnDirCIs{Index}(j,:) = ones(1,2)*NaN;
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
    
%     figure(DirUnDirComparisonFig);
%     subplot(3, ceil(length(DirBout_Stats)/3),Index);
%     hold on;
%     % plot(Distances, DistanceDirMeans(Index,:), 'rs-', 'LineWidth', 1);
%     errorbar(Distances, DistanceDirMeans(Index,:), DistanceDirSEMs(Index,:), 'rs-', 'LineWidth', 1);
%     % plot(repmat(Distances(:), 1, 2)', DistanceDirCIs{Index}', 'r')
%     
%     errorbar(Distances, DistanceUnDirMeans(Index,:), DistanceUnDirSEMs(Index,:), 'bs-', 'LineWidth', 1);
%     % plot(Distances, DistanceUnDirMeans(Index,:), 'bs-', 'LineWidth', 1);
%     % plot(repmat(Distances(:), 1, 2)', DistanceUnDirCIs{Index}', 'b');
%     plot(repmat(Distances(:),1,2), repmat(DistanceUnDirCIs{Index}(end,:), length(Distances), 1), 'k--');
%     
%     xlabel('Distance from female (cm)');
%     axis tight;
%     Temp = axis;
%     Temp = [-5 250 0.95*Temp(3) 1.05*Temp(4)];
%     axis(Temp);
%     ylabel(LabelString);
%     title({[BirdNames{i}, ': r=', num2str(Distance_DirSong_Corr(i,1)), '; p=', num2str(Distance_DirSong_Corr(i,2)), ';p=', num2str(KWP_DirSong{i})]; [BirdNames{i}, ': r=', num2str(Distance_UnDirSong_Corr(i,1)), '; p=', num2str(Distance_UnDirSong_Corr(i,2)), ';p=', num2str(KWP_UnDirSong{i})]}, 'Color', 'k', 'FontSize', 8, 'FontWeight', 'normal');
    % legend('D', 'DUN', 'UN');
end

OutputStats.DistanceDirMeans = DistanceDirMeans;
OutputStats.DistanceUnDirMeans = DistanceUnDirMeans;

OutputStats.DistanceDirSEMs = DistanceDirSEMs;
OutputStats.DistanceUnDirSEMs = DistanceUnDirSEMs;

%OutputStats.DistanceDirCIs = DistanceDirCIs;
OutputStats.DistanceUnDirCIs = DistanceUnDirCIs;

OutputStats.Distance_DirSong_Corr = Distance_DirSong_Corr;
OutputStats.Distance_UnDirSong_Corr = Distance_UnDirSong_Corr;

OutputStats.KWP_DirSong = KWP_DirSong;
OutputStats.KWP_UnDirSong = KWP_UnDirSong;

OutputStats.L0Dir_Undir_PValue = L0Dir_Undir_PValue;

OutputStats.LMFit_Intercept = LMFit_Intercept;
OutputStats.LMFit_DistanceContext_Interactions = LMFit_DistanceContext_Interactions;
OutputStats.LMFit_Distance = LMFit_Distance;
OutputStats.LMFit_Context = LMFit_Context;

OutputStats.LMFit_ModelPValue = LMFit_ModelPValue;

OutputStats.LMFit_ZScore_Intercept = LMFit_ZScore_Intercept;
OutputStats.LMFit_ZScore_DistanceContext_Interactions = LMFit_ZScore_DistanceContext_Interactions;
OutputStats.LMFit_ZScore_Distance = LMFit_ZScore_Distance;
OutputStats.LMFit_ZScore_Context = LMFit_ZScore_Context;

OutputStats.GroupData = GroupData;
OutputStats.Models = Models;

disp('Done with stats and plots');
