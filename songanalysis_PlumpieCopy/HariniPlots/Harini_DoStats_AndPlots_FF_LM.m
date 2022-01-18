function [OutputStats] = Harini_DoStats_AndPlots_FF_LM(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNot, SyllIndex)

% Now to do the stats
% 1. Compare L0 directed and undirected
% 2. For the ones where L0 directed and undirected are significantly
% different, have to do a linear correlation and see if there is a
% correlation or not
% 3. Plot individual L0 directed and undirected comparisons
% 4. Plot average % change across all birds over distance

Colours = 'rgbcmk';
CV = @(x) std(x)/mean(x); % a function handle for calculating the CV

for i = 1:length(DirBout_Stats)

    disp(['Bird #', num2str(i), ': ', BirdNames{i}]);
    if (~isfield(UndirBout_Stats{i}{end}, [StructString, '_SyllLabel']))
        disp('No syllables');
        continue;
    end
    for SyllIndex = 1:length(eval(['UndirBout_Stats{i}{end}.', StructString, '_SyllLabel'])),
        GroupData{i}{SyllIndex} = [];
        GroupData_CV{i}{SyllIndex} = [];
        LMGroupData{i}{SyllIndex} = [];
    end
    % GroupData has 3 columns, the first one has the data corresponding to
    % the variable that that I'm looking at like TotalINNumber, etc. 2nd
    % column has the index of the distance at which the recordings were
    % done (0 for L0, 1 for L1 ... 6 for UN) and 3rd column has 1 for Dir,
    % 2 for Dir/Undir - DUN - and 3 for Undir
    
    for j = 1:length(DirBout_Stats{i}),
        if (~isempty(DirBout_Stats{i}{j}))
            for SyllIndex = 1:length(eval(['DirBout_Stats{i}{j}.', StructString])),
                if (eval(['length(~isnan(DirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)))']) >= MinTrialNo)
                    GroupData{i}{SyllIndex} = [GroupData{i}{SyllIndex}; [eval(['DirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)']) ones(size(eval(['DirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)'])))*j ones(size(eval(['DirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)'])))*1]];
                    TempFF = eval(['DirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)']);
                    TempFF_CV = jackknife(CV, TempFF);
                    % GroupData_CV{i}{SyllIndex} = [GroupData_CV{i}{SyllIndex}; [std(TempFF)/mean(TempFF) ones(size(eval(['DirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)'])))*j ones(size(eval(['DirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)'])))*1]];
                    % LMGroupData{i}{SyllIndex} = [LMGroupData{i}{SyllIndex}; [TempFF_CV(:) ones(size(eval(['DirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)'])))*j ones(size(eval(['DirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)'])))*1]];
                    LMGroupData{i}{SyllIndex} = [LMGroupData{i}{SyllIndex}; [std(TempFF)/mean(TempFF) j 1]];
                else
                    LMGroupData{i}{SyllIndex} = [LMGroupData{i}{SyllIndex}; [ones(MinTrialNo,1)*NaN ones(MinTrialNo,1)*j ones(MinTrialNo,1)*1]];
                end
            end
        else
            LMGroupData{i}{SyllIndex} = [LMGroupData{i}{SyllIndex}; [ones(MinTrialNo,1)*NaN ones(MinTrialNo,1)*j ones(MinTrialNo,1)*1]];
        end
    end
    
    for j = 1:length(DirUndirBout_Stats{i}),
        if (~isempty(DirUndirBout_Stats{i}{j}))
            for SyllIndex = 1:length(eval(['DirUndirBout_Stats{i}{j}.', StructString])),
                if (eval(['length(~isnan(DirUndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)))']) >= MinTrialNo)
                    GroupData{i}{SyllIndex} = [GroupData{i}{SyllIndex}; [eval(['DirUndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)']) ones(size(eval(['DirUndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)'])))*j ones(size(eval(['DirUndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)'])))*2]];
                    TempFF = eval(['DirUndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)']);
                    TempFF_CV = jackknife(CV, TempFF);
                    GroupData_CV{i}{SyllIndex} = [GroupData_CV{i}{SyllIndex}; [TempFF_CV(:) ones(size(eval(['DirUndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)'])))*j ones(size(eval(['DirUndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)'])))*2]];
                    % LMGroupData{i}{SyllIndex} = [LMGroupData{i}{SyllIndex}; [TempFF_CV(:) ones(size(eval(['DirUndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)'])))*j ones(size(eval(['DirUndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)'])))*2]];
                    LMGroupData{i}{SyllIndex} = [LMGroupData{i}{SyllIndex}; [std(TempFF)/mean(TempFF) j 2]];
                else
                    LMGroupData{i}{SyllIndex} = [LMGroupData{i}{SyllIndex}; [ones(MinTrialNo,1)*NaN ones(MinTrialNo,1)*j ones(MinTrialNo,1)*2]];
                end
            end
        else
            LMGroupData{i}{SyllIndex} = [LMGroupData{i}{SyllIndex}; [ones(MinTrialNo,1)*NaN ones(MinTrialNo,1)*j ones(MinTrialNo,1)*2]];
        end
    end
    
    for j = 1:length(UndirBout_Stats{i}),
        if (~isempty(UndirBout_Stats{i}{j}))
            for SyllIndex = 1:length(eval(['UndirBout_Stats{i}{j}.', StructString])),
                if (eval(['length(~isnan(UndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)))']) >= MinTrialNo)
                    GroupData{i}{SyllIndex} = [GroupData{i}{SyllIndex}; [eval(['UndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)']) ones(size(eval(['UndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)'])))*j ones(size(eval(['UndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)'])))*3]];
                    TempFF = eval(['UndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)']);
                    TempFF_CV = jackknife(CV, TempFF);
                    GroupData_CV{i}{SyllIndex} = [GroupData_CV{i}{SyllIndex}; [TempFF_CV(:) ones(size(eval(['UndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)'])))*j ones(size(eval(['UndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)'])))*3]];
                    % LMGroupData{i}{SyllIndex} = [LMGroupData{i}{SyllIndex}; [TempFF_CV(:) ones(size(eval(['UndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)'])))*j ones(size(eval(['UndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)'])))*3]];
                    LMGroupData{i}{SyllIndex} = [LMGroupData{i}{SyllIndex}; [std(TempFF)/mean(TempFF) j 3]];
                else
                    LMGroupData{i}{SyllIndex} = [LMGroupData{i}{SyllIndex}; [ones(MinTrialNo,1)*NaN ones(MinTrialNo,1)*j ones(MinTrialNo,1)*3]];
                end
            end
        else
            LMGroupData{i}{SyllIndex} = [LMGroupData{i}{SyllIndex}; [ones(MinTrialNo,1)*NaN ones(MinTrialNo,1)*j ones(MinTrialNo,1)*3]];
        end
    end

    % Remove NaN values
    for SyllIndex = 1:length(GroupData{i}),
        NaNIndices = find(isnan(GroupData{i}{SyllIndex}(:,1)));
        if (~isempty(NaNIndices))
            GroupData{i}{SyllIndex}(NaNIndices,:) = [];
        end
        
        NaNIndices = find(isnan(GroupData_CV{i}{SyllIndex}(:,1)));
        if (~isempty(NaNIndices))
            GroupData_CV{i}{SyllIndex}(NaNIndices,:) = [];
        end
            
        % Run correlations - pearsons on dir data and undir data
    
        AllDirSongIndices = find(GroupData{i}{SyllIndex}(:,3) == 1);
        [r, p] = corrcoef(GroupData{i}{SyllIndex}(AllDirSongIndices, 1), GroupData{i}{SyllIndex}(AllDirSongIndices, 2));
        Distance_DirSong_Corr{i}(SyllIndex,:) = [r(1,2) p(1,2)];
        
        if (ParametricOrNot == 0)
            [KWP_DirSong{i}(SyllIndex), anovatab, KWStats_DirSong{i}{SyllIndex}] = kruskalwallis(GroupData{i}{SyllIndex}(AllDirSongIndices, 1), GroupData{i}{SyllIndex}(AllDirSongIndices, 2), 'off');
        else
            [KWP_DirSong{i}(SyllIndex), anovatab, KWStats_DirSong{i}{SyllIndex}] = anova1(GroupData{i}{SyllIndex}(AllDirSongIndices, 1), GroupData{i}{SyllIndex}(AllDirSongIndices, 2), 'off');
            [VarTestP_DirSong{i}(SyllIndex), VarTestStats_DirSong{i}{SyllIndex}] = vartestn(GroupData{i}{SyllIndex}(AllDirSongIndices, 1), GroupData{i}{SyllIndex}(AllDirSongIndices, 2), 'display', 'off');
        end        

        AllUnDirSongIndices = find((GroupData{i}{SyllIndex}(:,3) == 3) & (GroupData{i}{SyllIndex}(:,2) < 6));
        [r, p] = corrcoef(GroupData{i}{SyllIndex}(AllUnDirSongIndices, 1), GroupData{i}{SyllIndex}(AllUnDirSongIndices, 2));
        Distance_UnDirSong_Corr{i}(SyllIndex,:) = [r(1,2) p(1,2)];
        if (~isempty(AllUnDirSongIndices))
            if (ParametricOrNot == 0)
                [KWP_UnDirSong{i}(SyllIndex), anovatab, KWStats_UnDirSong{i}{SyllIndex}] = kruskalwallis(GroupData{i}{SyllIndex}(AllUnDirSongIndices, 1), GroupData{i}{SyllIndex}(AllUnDirSongIndices, 2), 'off');
            else
                [KWP_UnDirSong{i}(SyllIndex), anovatab, KWStats_UnDirSong{i}{SyllIndex}] = anova1(GroupData{i}{SyllIndex}(AllUnDirSongIndices, 1), GroupData{i}{SyllIndex}(AllUnDirSongIndices, 2), 'off');
                [VarTestP_UnDirSong{i}(SyllIndex), VarTestStats_UnDirSong{i}{SyllIndex}] = vartestn(GroupData{i}{SyllIndex}(AllUnDirSongIndices, 1), GroupData{i}{SyllIndex}(AllUnDirSongIndices, 2), 'display', 'off');
            end
        else
            KWP_UnDirSong{i}(SyllIndex) = NaN;
            VarTestP_UnDirSong{i}(SyllIndex) = NaN;
        end
        L0Dir_Indices = find((GroupData{i}{SyllIndex}(:,2) == 1) & (GroupData{i}{SyllIndex}(:,3) == 1));
        UnDir_Indices = find((GroupData{i}{SyllIndex}(:,2) == 6) & (GroupData{i}{SyllIndex}(:,3) == 3));

        if ((length(L0Dir_Indices) >= MinTrialNo) && (length(UnDir_Indices) >= MinTrialNo))
            switch ParametricOrNot
                case 0
                    L0Dir_Undir_PValue{i}(SyllIndex) = kruskalwallis([GroupData{i}{SyllIndex}(L0Dir_Indices,1); GroupData{i}{SyllIndex}(UnDir_Indices,1)], [ones(size(L0Dir_Indices(:)))*1; ones(size(UnDir_Indices(:)))*2], 'off');
                    % L0Dir_Undir_PValue(i) = signrank(GroupData{i}(L0Dir_Indices,1), GroupData{i}(UnDir_Indices,1));
                case 1
                    [HypothesisValue, L0Dir_Undir_PValue{i}(SyllIndex)] = ttest2(GroupData{i}{SyllIndex}(L0Dir_Indices,1), GroupData{i}{SyllIndex}(UnDir_Indices,1));
                    [HypothesisValue, L0Dir_Undir_VarTestPValue{i}(SyllIndex)] = vartest2(GroupData{i}{SyllIndex}(L0Dir_Indices,1), GroupData{i}{SyllIndex}(UnDir_Indices,1));
            end
        else
            L0Dir_Undir_PValue{i}(SyllIndex) = NaN;
            L0Dir_Undir_VarTestPValue{i}(SyllIndex) = NaN;
        end
    
        DisplayString = ['P value for ', StructString, ': '];
        if (ParametricOrNot == 0)
            DisplayString = [DisplayString, 'Kruskalwallis test '];
        else
            DisplayString = [DisplayString, 't-test '];
        end
    
        DisplayString = [DisplayString, 'comparing mean FF of L0 Directed songs and Undirected songs for syll ', eval(['UndirBout_Stats{i}{end}.', StructString, '_SyllLabel(', num2str(SyllIndex), ')']), ' = ', num2str(L0Dir_Undir_PValue{i}(SyllIndex))];
        DisplayString = [DisplayString, 'comparing var FF of L0 Directed songs and Undirected songs for syll ', eval(['UndirBout_Stats{i}{end}.', StructString, '_SyllLabel(', num2str(SyllIndex), ')']), ' = ', num2str(L0Dir_Undir_VarTestPValue{i}(SyllIndex))];
        
        if (mean(GroupData{i}{SyllIndex}(L0Dir_Indices,1)) > mean(GroupData{i}{SyllIndex}(UnDir_Indices,1)))
            DisplayString = [DisplayString, '; Mean FF: Dir > Undir'];
        else
            if (mean(GroupData{i}{SyllIndex}(L0Dir_Indices,1)) < mean(GroupData{i}{SyllIndex}(UnDir_Indices,1)))
                DisplayString = [DisplayString, '; MeanFF: Dir < Undir'];
            else
                DisplayString = [DisplayString, '; MeanFF: Dir = Undir'];
            end
        end
        
        if ((std(GroupData{i}{SyllIndex}(L0Dir_Indices,1))/mean(GroupData{i}{SyllIndex}(L0Dir_Indices,1))) > (std(GroupData{i}{SyllIndex}(UnDir_Indices,1))/mean(GroupData{i}{SyllIndex}(UnDir_Indices,1))))
            DisplayString = [DisplayString, '; VarFF: Dir > Undir'];
        else
            if ((std(GroupData{i}{SyllIndex}(L0Dir_Indices,1))/mean(GroupData{i}{SyllIndex}(L0Dir_Indices,1))) < (std(GroupData{i}{SyllIndex}(UnDir_Indices,1))/mean(GroupData{i}{SyllIndex}(UnDir_Indices,1))))
                DisplayString = [DisplayString, '; VarFF: Dir < Undir'];
            else
                DisplayString = [DisplayString, '; VarFF: Dir = Undir'];
            end
        end
        
        disp(DisplayString);

        if ((~isempty(L0Dir_Indices)) && (~isempty(UnDir_Indices)))
            L0Dir_MeanValues{i}(SyllIndex) = mean(GroupData{i}{SyllIndex}(L0Dir_Indices,1));
            L0Dir_MedianValues{i}(SyllIndex) = median(GroupData{i}{SyllIndex}(L0Dir_Indices,1));
            L0Dir_MaxValues{i}(SyllIndex) = max(GroupData{i}{SyllIndex}(L0Dir_Indices,1));
            UnDir_MeanValues{i}(SyllIndex) = mean(GroupData{i}{SyllIndex}(UnDir_Indices,1));
            UnDir_MedianValues{i}(SyllIndex) = median(GroupData{i}{SyllIndex}(UnDir_Indices,1));
            UnDir_MaxValues{i}(SyllIndex) = max(GroupData{i}{SyllIndex}(UnDir_Indices,1));
            L0Dir_CVValues{i}(SyllIndex) = std(GroupData{i}{SyllIndex}(L0Dir_Indices,1))/mean(GroupData{i}{SyllIndex}(L0Dir_Indices,1));
            UnDir_CVValues{i}(SyllIndex) = std(GroupData{i}{SyllIndex}(UnDir_Indices,1))/mean(GroupData{i}{SyllIndex}(UnDir_Indices,1));
        else
            L0Dir_MeanValues{i}(SyllIndex) = NaN;
            L0Dir_MedianValues{i}(SyllIndex) = NaN;
            L0Dir_MaxValues{i}(SyllIndex) = NaN;
            UnDir_MeanValues{i}(SyllIndex) = NaN;
            UnDir_MedianValues{i}(SyllIndex) = NaN;
            UnDir_MaxValues{i}(SyllIndex) = NaN;
            L0Dir_CVValues{i}(SyllIndex) = NaN;
            UnDir_CVValues{i}(SyllIndex) = NaN;
        end
        
        TempDataIndices = find(LMGroupData{i}{SyllIndex}(:,3) == 2);
        Temp2DataIndices = find(LMGroupData{i}{SyllIndex}(:,2) == 6);
        LMModelDataIndices = setdiff(1:1:size(LMGroupData{i}{SyllIndex},1), union(TempDataIndices, Temp2DataIndices));
        
        ContextValues = {'D' 'DUN' 'UN'};
        ContextNumericalValues = [0 2 1];

        FF_CV = LMGroupData{i}{SyllIndex}(LMModelDataIndices,1);
        Contexts = ContextValues(LMGroupData{i}{SyllIndex}(LMModelDataIndices,3));
        DistanceValues = Distances(LMGroupData{i}{SyllIndex}(LMModelDataIndices,2));

        Context_NumericalValues = ContextNumericalValues(LMGroupData{i}{SyllIndex}(LMModelDataIndices, 3));

        FF_CV = FF_CV(:);
        Contexts = Contexts(:);
        DistanceValues = DistanceValues(:);

        LM_DataTable = table(DistanceValues, Contexts, FF_CV);
        mdl = fitlm(LM_DataTable, 'interactions', 'CategoricalVars', [2]);

        mdl3 = fitlm(([zscore(DistanceValues(:)) Context_NumericalValues(:)]), (FF_CV(:) - nanmean(FF_CV))/nanstd(FF_CV), 'interactions', 'CategoricalVars', [2]);
        Temp_Estimates_PValues = mdl.Coefficients{:, {'Estimate', 'pValue'}};

        if (size(Temp_Estimates_PValues,1) < 4)
            LMFit_ModelPValue{i}(SyllIndex) = NaN;
        else
            Temp_Model_Summary = anova(mdl, 'summary');
            LMFit_ModelPValue{i}(SyllIndex) = table2array(Temp_Model_Summary(2,5));
        end
        
        if (size(Temp_Estimates_PValues,1) == 4)
            LMFit_Intercept{i}(SyllIndex,:) = Temp_Estimates_PValues(1,:);
            LMFit_Distance{i}(SyllIndex,:) = Temp_Estimates_PValues(2,:);
            LMFit_Context{i}(SyllIndex,:) = Temp_Estimates_PValues(3,:);
            LMFit_DistanceContext_Interactions{i}(SyllIndex,:) = Temp_Estimates_PValues(4,:);
        else
            if (size(Temp_Estimates_PValues,1) == 3)
                LMFit_Intercept{i}(SyllIndex,:) = Temp_Estimates_PValues(1,:);
                LMFit_Distance{i}(SyllIndex,:) = Temp_Estimates_PValues(2,:);
                LMFit_Context{i}(SyllIndex,:) = Temp_Estimates_PValues(3,:);
                LMFit_DistanceContext_Interactions{i}(SyllIndex,:) = ones(1,2)*NaN;
            else
                LMFit_Intercept{i}(SyllIndex,:) = ones(1,2)*NaN;
                LMFit_Distance{i}(SyllIndex,:) = ones(1,2)*NaN;
                LMFit_Context{i}(SyllIndex,:) = ones(1,2)*NaN;
                LMFit_DistanceContext_Interactions{i}(SyllIndex,:) = ones(1,2)*NaN;
            end
        end
        
            
        Temp_Estimates_PValues = mdl3.Coefficients{:, {'Estimate', 'pValue'}};
        if (size(Temp_Estimates_PValues,1) == 4)
            LMFit_ZScore_Intercept{i}(SyllIndex,:) = Temp_Estimates_PValues(1,:);
            LMFit_ZScore_Distance{i}(SyllIndex,:) = Temp_Estimates_PValues(2,:);
            LMFit_ZScore_Context{i}(SyllIndex,:) = Temp_Estimates_PValues(3,:);
            LMFit_ZScore_DistanceContext_Interactions{i}(SyllIndex,:) = Temp_Estimates_PValues(4,:);
        else
            if (size(Temp_Estimates_PValues,1) == 3)
                LMFit_ZScore_Intercept{i}(SyllIndex,:) = Temp_Estimates_PValues(1,:);
                LMFit_ZScore_Distance{i}(SyllIndex,:) = Temp_Estimates_PValues(2,:);
                LMFit_ZScore_Context{i}(SyllIndex,:) = Temp_Estimates_PValues(3,:);
                LMFit_ZScore_DistanceContext_Interactions{i}(SyllIndex,:) = ones(1,2)*NaN;
            else
                LMFit_ZScore_Intercept{i}(SyllIndex,:) = ones(1,2)*NaN;
                LMFit_ZScore_Distance{i}(SyllIndex,:) = ones(1,2)*NaN;
                LMFit_ZScore_Context{i}(SyllIndex,:) = ones(1,2)*NaN;
                LMFit_ZScore_DistanceContext_Interactions{i}(SyllIndex,:) = ones(1,2)*NaN;
            end
        end

        Models{i}{SyllIndex}.mdl = mdl;
        Models{i}{SyllIndex}.ZScoremdl = mdl3;
        disp(mdl);
        disp(mdl3);
    end
end

% disp('Distance directed song correlations: ');
% disp(Distance_DirSong_Corr);
% 
% disp('Distance undirected song correlations: ');
% disp(Distance_UnDirSong_Corr);

disp('Group data comparisons');

figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [131 558 700 400]);
p = panel();
p.pack('h', {1/2 1/2});
p.de.margin = 10;

p(1).select();
hold on;
AllL0Dir_UnDirCVs = [cell2mat(L0Dir_CVValues); cell2mat(UnDir_CVValues); cell2mat(L0Dir_Undir_VarTestPValue)];
for i = 1:2,
    L0DirUnDirBar(i) = bar(i, nanmean(AllL0Dir_UnDirCVs(i,:)));
    set(L0DirUnDirBar(i), 'FaceColor', 'none', 'LineWidth', 1);
    errorbar(i, nanmean(AllL0Dir_UnDirCVs(i,:)), (nanstd(AllL0Dir_UnDirCVs(i,:))/sqrt(length(find(~isnan(AllL0Dir_UnDirCVs(i,:)))))), 'k', 'LineWidth', 1);
end

SigCVs = find(AllL0Dir_UnDirCVs(3,:) < 0.05);
plot(repmat([1.1; 1.9], 1, length(SigCVs)), AllL0Dir_UnDirCVs(1:2,SigCVs), 'ko-', 'MarkerFaceColor', 'k');

NonSigCVs = find(AllL0Dir_UnDirCVs(3,:) >= 0.05);
plot(repmat([1.1; 1.9], 1, length(NonSigCVs)), AllL0Dir_UnDirCVs(1:2,NonSigCVs), 'ko-');

set(gca, 'XTick', [1 2], 'XTickLabel', {'L0 Dir' 'Undir'});
axis auto;
ylabel(['CV ', LabelString]);

Group_L0Dir_UnDir_PValue.CV = signrank(AllL0Dir_UnDirCVs(1, find(~isnan(AllL0Dir_UnDirCVs(1,:)))), AllL0Dir_UnDirCVs(2, find(~isnan(AllL0Dir_UnDirCVs(1,:)))));
sigstar({[1 2]}, Group_L0Dir_UnDir_PValue.CV);
if (nanmean(AllL0Dir_UnDirCVs(1,:)) > nanmean(AllL0Dir_UnDirCVs(2,:)))
    disp(['Sign rank test comparing group mean CVs for L0 directed and Undirected, p = ', num2str(Group_L0Dir_UnDir_PValue.CV), '; Dir > Undir']);
else
    if (nanmean(AllL0Dir_UnDirCVs(1,:)) < nanmean(AllL0Dir_UnDirCVs(2,:)))
        disp(['Sign rank test comparing group mean CVs for L0 directed and Undirected, p = ', num2str(Group_L0Dir_UnDir_PValue.CV), '; Dir < Undir']);
    else
        disp(['Sign rank test comparing group mean CVs for L0 directed and Undirected, p = ', num2str(Group_L0Dir_UnDir_PValue.CV), '; Dir = Undir']);
    end
end

p(2).select();
hold on;
AllL0Dir_UnDirMeans = [cell2mat(L0Dir_MeanValues); cell2mat(UnDir_MeanValues); cell2mat(L0Dir_Undir_PValue)];
for i = 1:2,
    L0DirUnDirBar(i) = bar(i, nanmean(AllL0Dir_UnDirMeans(i,:)));
    set(L0DirUnDirBar(i), 'FaceColor', 'none', 'LineWidth', 1);
    errorbar(i, nanmean(AllL0Dir_UnDirMeans(i,:)), (nanstd(AllL0Dir_UnDirMeans(i,:))/sqrt(length(find(~isnan(AllL0Dir_UnDirMeans(i,:)))))), 'k', 'LineWidth', 1);
end

SigCVs = find(AllL0Dir_UnDirMeans(3,:) < 0.05);
plot(repmat([1.1; 1.9], 1, length(SigCVs)), AllL0Dir_UnDirMeans(1:2,SigCVs), 'ko-', 'MarkerFaceColor', 'k');

NonSigCVs = find(AllL0Dir_UnDirMeans(3,:) >= 0.05);
plot(repmat([1.1; 1.9], 1, length(NonSigCVs)), AllL0Dir_UnDirMeans(1:2,NonSigCVs), 'ko-');

set(gca, 'XTick', [1 2], 'XTickLabel', {'L0 Dir' 'Undir'});
axis auto;
ylabel(['CV ', LabelString]);

Group_L0Dir_UnDir_PValue.Mean = signrank(AllL0Dir_UnDirMeans(1, find(~isnan(AllL0Dir_UnDirMeans(1,:)))), AllL0Dir_UnDirMeans(2, find(~isnan(AllL0Dir_UnDirMeans(1,:)))));
sigstar({[1 2]}, Group_L0Dir_UnDir_PValue.Mean);
if (nanmean(AllL0Dir_UnDirMeans(1,:)) > nanmean(AllL0Dir_UnDirMeans(2,:)))
    disp(['Sign rank test comparing group mean CVs for L0 directed and Undirected, p = ', num2str(Group_L0Dir_UnDir_PValue.Mean), '; Dir > Undir']);
else
    if (nanmean(AllL0Dir_UnDirMeans(1,:)) < nanmean(AllL0Dir_UnDirMeans(2,:)))
        disp(['Sign rank test comparing group means for L0 directed and Undirected, p = ', num2str(Group_L0Dir_UnDir_PValue.Mean), '; Dir < Undir']);
    else
        disp(['Sign rank test comparing group means for L0 directed and Undirected, p = ', num2str(Group_L0Dir_UnDir_PValue.Mean), '; Dir = Undir']);
    end
end

        

TotalNumSylls = length(cell2mat(KWP_DirSong));

DistanceDirMeans = [];
DistanceDirUnDirMeans = [];
DistanceUnDirMeans = [];
Symbols = 'so^+hp><*dxv.so';

Index = 0;

DirUnDirCVComparisonFig = figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [188 208 1600 800]);

DirUnDirMeanComparisonFig = figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [188 208 1600 800]);

Index = 0;


for i = 1:length(BirdNames),
    for SyllIndex = 1:length(GroupData{i}),
        Index = Index + 1;
        AllDirSongIndices = find(GroupData{i}{SyllIndex}(:,3) == 1);

        % First for directed songs
        for j = 1:6,
            DirSongIndices = find((GroupData{i}{SyllIndex}(:,3) == 1) & (GroupData{i}{SyllIndex}(:,2) == j));
            DirSong_NumTrials(i,j) = length(DirSongIndices);
            if (length(DirSongIndices) >= MinTrialNo)
                DistanceDirMeans(Index, j) = mean(GroupData_CV{i}{SyllIndex}(DirSongIndices,1));
                DistanceDirSEMs(Index, j) = std(GroupData_CV{i}{SyllIndex}(DirSongIndices,1))/sqrt(length(DirSongIndices));
                DistanceDirCIs{Index}(j,:) = bootci(10000,@mean, GroupData{i}{SyllIndex}(DirSongIndices,1));
                DistanceDirCVs(Index,j) = std(GroupData{i}{SyllIndex}(DirSongIndices,1))/mean(GroupData{i}{SyllIndex}(DirSongIndices,1));
                DistanceDirCVCIs{Index}(j,:) = bootci(10000, CV, GroupData{i}{SyllIndex}(DirSongIndices, 1));
            else
                DistanceDirMeans(Index, j) = NaN;
                DistanceDirSEMs(Index,j) = NaN;
                DistanceDirCIs{Index}(j,:) = ones(1,2)*NaN;
                DistanceDirCVs(Index, j) = NaN;
                DistanceDirCVCIs{Index}(j,:) = ones(1,2)*NaN;
            end

            DirUnDirSongIndices = find((GroupData{i}{SyllIndex}(:,3) == 2) & (GroupData{i}{SyllIndex}(:,2) == j));
            DirUnDirSong_NumTrials(i,j) = length(DirUnDirSongIndices);
            if (length(DirUnDirSongIndices) >= MinTrialNo)
                DistanceDirUnDirMeans(Index, j) = mean(GroupData{i}{SyllIndex}(DirUnDirSongIndices,1));
                DistanceDirUnDirSEMs(Index, j) = std(GroupData{i}{SyllIndex}(DirUnDirSongIndices,1))/sqrt(length(DirUnDirSongIndices));
                DistanceDirUnDirCIs{Index}(j,:) = bootci(10000, @mean, GroupData{i}{SyllIndex}(DirUnDirSongIndices,1));
                DistanceDirUnDirCVs(Index,j) = std(GroupData{i}{SyllIndex}(DirUnDirSongIndices,1))/mean(GroupData{i}{SyllIndex}(DirUnDirSongIndices,1));
                DistanceDirUnDirCVCIs{Index}(j,:) = bootci(10000, CV, GroupData{i}{SyllIndex}(DirUnDirSongIndices, 1));
            else
                DistanceDirUnDirMeans(Index, j) = NaN;
                DistanceDirUnDirSEMs(Index, j) = NaN;
                DistanceDirUnDirCIs{Index}(j,:) = ones(1,2)*NaN;
                DistanceDirUnDirCVs(Index, j) = NaN;
                DistanceDirUnDirCVCIs{Index}(j,:) = ones(1,2)*NaN;
            end

            UnDirSongIndices = find((GroupData{i}{SyllIndex}(:,3) == 3) & (GroupData{i}{SyllIndex}(:,2) == j));
            UnDirSong_NumTrials(i,j) = length(UnDirSongIndices);
            if (length(UnDirSongIndices) >= MinTrialNo)
                DistanceUnDirMeans(Index, j) = mean(GroupData{i}{SyllIndex}(UnDirSongIndices,1));
                DistanceUnDirSEMs(Index, j) = std(GroupData{i}{SyllIndex}(UnDirSongIndices,1))/sqrt(length(UnDirSongIndices));
                DistanceUnDirCIs{Index}(j,:) = bootci(10000, @mean, GroupData{i}{SyllIndex}(UnDirSongIndices,1));
                DistanceUnDirCVs(Index,j) = std(GroupData{i}{SyllIndex}(UnDirSongIndices,1))/mean(GroupData{i}{SyllIndex}(UnDirSongIndices,1));
                DistanceUnDirCVCIs{Index}(j,:) = bootci(10000, CV, GroupData{i}{SyllIndex}(UnDirSongIndices, 1));
            else
                DistanceUnDirMeans(Index, j) = NaN;
                DistanceUnDirSEMs(Index, j) = NaN;
                DistanceUnDirCIs{Index}(j,:) = ones(1,2)*NaN;
                DistanceUnDirCVs(Index, j) = NaN;
                DistanceUnDirCVCIs{Index}(j,:) = ones(1,2)*NaN;
            end
        end

        figure(DirUnDirCVComparisonFig);
        subplot(ceil(TotalNumSylls/6), 6, Index);
        hold on;
        plot(Distances, DistanceDirCVs(Index,:), 'rs-', 'LineWidth', 1);
        plot(repmat(Distances(:), 1, 2)', DistanceDirCVCIs{Index}', 'r')

    %     plot(Distances, DistanceDirUnDirMeans(Index,:), 'ks-', 'LineWidth', 1);
    %     plot(repmat(Distances(:), 1, 2)', DistanceDirUnDirCIs{Index}', 'k')
    %     
        plot(Distances, DistanceUnDirCVs(Index,:), 'bs-', 'LineWidth', 1);
        plot(repmat(Distances(:), 1, 2)', DistanceUnDirCVCIs{Index}', 'b');
        plot(repmat(Distances(:),1,2), repmat(DistanceUnDirCVCIs{Index}(end,:), length(Distances), 1), 'k--');

                
        figure(DirUnDirMeanComparisonFig);
        subplot(ceil(TotalNumSylls/6), 6, Index);
        hold on;
        plot(Distances, DistanceDirMeans(Index,:), 'rs-', 'LineWidth', 1);
        plot(repmat(Distances(:), 1, 2)', DistanceDirCIs{Index}', 'r')

    %     plot(Distances, DistanceDirUnDirMeans(Index,:), 'ks-', 'LineWidth', 1);
    %     plot(repmat(Distances(:), 1, 2)', DistanceDirUnDirCIs{Index}', 'k')
    %     
        plot(Distances, DistanceUnDirMeans(Index,:), 'bs-', 'LineWidth', 1);
        plot(repmat(Distances(:), 1, 2)', DistanceUnDirCIs{Index}', 'b');
        plot(repmat(Distances(:),1,2), repmat(DistanceUnDirCIs{Index}(end,:), length(Distances), 1), 'k--');
        
        %xlabel('Distance from female (cm)');
        
        axis tight;
        Temp = axis;
        Temp = [-5 250 0.95*Temp(3) 1.05*Temp(4)];
        axis(Temp);
        %ylabel(LabelString);
        title({[BirdNames{i}, ': r=', num2str(Distance_DirSong_Corr{i}(SyllIndex,1)), '; p=', num2str(Distance_DirSong_Corr{i}(SyllIndex,2)), '; p=', num2str(KWP_DirSong{i}(SyllIndex))]; [BirdNames{i}, ': r=', num2str(Distance_UnDirSong_Corr{i}(SyllIndex,1)), '; p=', num2str(Distance_UnDirSong_Corr{i}(SyllIndex,2)), '; p=', num2str(KWP_UnDirSong{i}(SyllIndex))]}, 'Color', 'k', 'FontSize', 8, 'FontWeight', 'normal');
        
        figure(DirUnDirCVComparisonFig);
        axis tight;
        Temp = axis;
        Temp = [-5 250 0.95*Temp(3) 1.05*Temp(4)];
        axis(Temp);
        %ylabel(LabelString);
        title({[BirdNames{i}, ': p=', num2str(VarTestP_DirSong{i}(SyllIndex))]; [BirdNames{i}, ': p=', num2str(VarTestP_UnDirSong{i}(SyllIndex))]}, 'Color', 'k', 'FontSize', 8, 'FontWeight', 'normal');

        % Set the value for directed songs and dir/undir songs for the
        % undirected to be equal to the undirected song means

        DistanceDirMeans(Index,end) = DistanceUnDirMeans(Index,end);
        DistanceDirUnDirMeans(Index,end) = DistanceUnDirMeans(Index,end);


        DistanceDirMeans(Index,:) = 100*DistanceDirMeans(Index,:)/DistanceDirMeans(Index,end);
        DistanceDirUnDirMeans(Index,:) = 100*DistanceDirUnDirMeans(Index,:)/DistanceDirUnDirMeans(Index,end);
        DistanceUnDirMeans(Index,:) = 100*DistanceUnDirMeans(Index,:)/DistanceUnDirMeans(Index,end);

%         if (~isempty(find(BirdsWithSigDiffs == i)))
%             figure(GroupDataFig);
%             if (L0Dir_MeanValues(i) > UnDir_MeanValues(i))
%                 p(1).select();
%                 plot(Distances, DistanceDirMeans(Index,:), ['r', Symbols(i), '-']);
% 
%                 p(2).select();
%                 plot(Distances, DistanceDirUnDirMeans(Index,:), ['r', Symbols(i), '-']);
% 
%                 p(3).select();
%                 plot(Distances, DistanceUnDirMeans(Index,:), ['r', Symbols(i), '-']);
%             else
%                 if (L0Dir_MeanValues(i) < UnDir_MeanValues(i))
%                     p(1).select();
%                     plot(Distances, DistanceDirMeans(Index,:), ['b', Symbols(i), '-']);
% 
%                     p(2).select();
%                     plot(Distances, DistanceDirUnDirMeans(Index,:), ['b', Symbols(i), '-']);
% 
%                     p(3).select();
%                     plot(Distances, DistanceUnDirMeans(Index,:), ['b', Symbols(i), '-']);
%                 else
%                     p(1).select();
%                     title('Directed songs');
%                     plot(Distances, DistanceDirMeans(Index,:), ['k', Symbols(i), '-']);
% 
%                     p(2).select();
%                     title('DUN songs');
%                     plot(Distances, DistanceDirUnDirMeans(Index,:), ['k', Symbols(i), '-']);
% 
%                     p(3).select();
%                     title('Undir songs');
%                     plot(Distances, DistanceUnDirMeans(Index,:), ['k', Symbols(i), '-']);
%                 end
%             end
%         end
    end
end


% figure(GroupDataFig);
% for i = 1:3,
%     p(i).select();
%     if (i == 2)
%         ylabel(['% normalized change in ', LabelString]);
%     end
%     if (i == 3)
%         xlabel('Distance to female (cm)');
%     end
% end
cd 
OutputStats.DistanceDirMeans = DistanceDirMeans;
OutputStats.DistanceUnDirMeans = DistanceUnDirMeans;

OutputStats.DistanceDirCVs = DistanceDirCVs;
OutputStats.DistanceUnDirCVs = DistanceUnDirCVs;

OutputStats.Distance_DirSong_Corr = Distance_DirSong_Corr;
OutputStats.Distance_UnDirSong_Corr = Distance_UnDirSong_Corr;

OutputStats.KWP_DirSong = KWP_DirSong;
OutputStats.KWP_UnDirSong = KWP_UnDirSong;

OutputStats.L0Dir_Undir_PValue = L0Dir_Undir_PValue;

disp('Done with stats and plots');
