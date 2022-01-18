function [OutputStats] = Harini_DoStats_AndPlots_FF_LogAmp(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames, ParametricOrNot, SyllIndex)

% Now to do the stats
% 1. Compare L0 directed and undirected
% 2. For the ones where L0 directed and undirected are significantly
% different, have to do a linear correlation and see if there is a
% correlation or not
% 3. Plot individual L0 directed and undirected comparisons
% 4. Plot average % change across all birds over distance

Colours = 'rgbcmk';

for i = 1:length(DirBout_Stats)
    figure;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [188 208 400 800]);

    disp(['Bird #', num2str(i), ': ', BirdNames{i}]);
    if (~isfield(UndirBout_Stats{i}{end}, [StructString, '_SyllLabel']))
        disp('No syllables');
        continue;
    end
    for SyllIndex = 1:length(eval(['UndirBout_Stats{i}{end}.', StructString, '_SyllLabel'])),
        GroupData{i}{SyllIndex} = [];
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
                end
            end
        end
    end
    
    for j = 1:length(DirUndirBout_Stats{i}),
        if (~isempty(DirUndirBout_Stats{i}{j}))
            for SyllIndex = 1:length(eval(['DirUndirBout_Stats{i}{j}.', StructString])),
                if (eval(['length(~isnan(DirUndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)))']) >= MinTrialNo)
                    GroupData{i}{SyllIndex} = [GroupData{i}{SyllIndex}; [eval(['DirUndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)']) ones(size(eval(['DirUndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)'])))*j ones(size(eval(['DirUndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)'])))*2]];
                end
            end
        end
    end
    
    for j = 1:length(UndirBout_Stats{i}),
        if (~isempty(UndirBout_Stats{i}{j}))
            for SyllIndex = 1:length(eval(['UndirBout_Stats{i}{j}.', StructString])),
                if (eval(['length(~isnan(UndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)))']) >= MinTrialNo)
                    GroupData{i}{SyllIndex} = [GroupData{i}{SyllIndex}; [eval(['UndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)']) ones(size(eval(['UndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)'])))*j ones(size(eval(['UndirBout_Stats{i}{j}.', StructString, '{', num2str(SyllIndex), '}(:)'])))*3]];
                end
            end
        end
    end

    % Remove NaN values
    for SyllIndex = 1:length(GroupData{i}),
        NaNIndices = find(isnan(GroupData{i}{SyllIndex}(:,1)));
        if (~isempty(NaNIndices))
            GroupData{i}{SyllIndex}(NaNIndices,:) = [];
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
    
        switch StructString
            case 'LogAmplitude'
                DisplayString = [DisplayString, 'comparing L0 Directed songs and Undirected songs for syll ', eval(['UndirBout_Stats{i}{end}.', StructString, '_SyllLabel(', num2str(SyllIndex), ')']), ' = ', num2str(L0Dir_Undir_PValue{i}(SyllIndex))];
            case 'FF'
                DisplayString = [DisplayString, 'comparing L0 Directed songs and Undirected songs for syll ', eval(['UndirBout_Stats{i}{end}.', StructString, '_SyllLabel(', num2str(SyllIndex), ')']), ' = ', num2str(L0Dir_Undir_VarTestPValue{i}(SyllIndex))];
        end    
        
        switch StructString
            case 'LogAmplitude'
                if (mean(GroupData{i}{SyllIndex}(L0Dir_Indices,1)) > mean(GroupData{i}{SyllIndex}(UnDir_Indices,1)))
                    DisplayString = [DisplayString, '; Dir > Undir'];
                else
                    if (mean(GroupData{i}{SyllIndex}(L0Dir_Indices,1)) < mean(GroupData{i}{SyllIndex}(UnDir_Indices,1)))
                        DisplayString = [DisplayString, '; Dir < Undir'];
                    else
                        DisplayString = [DisplayString, '; Dir = Undir'];
                    end
                end
            case 'FF'
                if ((std(GroupData{i}{SyllIndex}(L0Dir_Indices,1))/mean(GroupData{i}{SyllIndex}(L0Dir_Indices,1))) > (std(GroupData{i}{SyllIndex}(UnDir_Indices,1))/mean(GroupData{i}{SyllIndex}(UnDir_Indices,1))))
                    DisplayString = [DisplayString, '; Dir > Undir'];
                else
                    if ((std(GroupData{i}{SyllIndex}(L0Dir_Indices,1))/mean(GroupData{i}{SyllIndex}(L0Dir_Indices,1))) < (std(GroupData{i}{SyllIndex}(UnDir_Indices,1))/mean(GroupData{i}{SyllIndex}(UnDir_Indices,1))))
                        DisplayString = [DisplayString, '; Dir < Undir'];
                    else
                        DisplayString = [DisplayString, '; Dir = Undir'];
                    end
                end
        end
        disp(DisplayString);

        % Now plot the distributions and means for all birds across all
        % distances for directed songs
        subplot(length(GroupData{i}), 1, SyllIndex);
        hold on;
        Edges = linspace(min(GroupData{i}{SyllIndex}(:,1)), max(GroupData{i}{SyllIndex}(:,1)), 5);

        for j = min(GroupData{i}{SyllIndex}(:,2)):max(GroupData{i}{SyllIndex}(:,2)),
            DistIndices = find((GroupData{i}{SyllIndex}(:,2) == j) & (GroupData{i}{SyllIndex}(:,3) == 1));
            plot(Edges, 100*histc(GroupData{i}{SyllIndex}(DistIndices,1), Edges)/length(DistIndices), [Colours(j),'o-']);
        end

        axis tight;
        Temp = axis;

        for j = min(GroupData{i}{SyllIndex}(:,2)):max(GroupData{i}{SyllIndex}(:,2)),
            DistIndices = find((GroupData{i}{SyllIndex}(:,2) == j) & (GroupData{i}{SyllIndex}(:,3) == 1));
            plot(mean(GroupData{i}{SyllIndex}(DistIndices,1)), (1+j*0.05)*1.05*Temp(4), [Colours(j), 's']);
            plot(mean(GroupData{i}{SyllIndex}(DistIndices,1)) + [std(GroupData{i}{SyllIndex}(DistIndices,1))/sqrt(length(DistIndices)) -std(GroupData{i}{SyllIndex}(DistIndices,1))/sqrt(length(DistIndices))], ones(1,2)*(1+j*0.05)*1.05*Temp(4), Colours(j));
        end

        axis tight;
        xlabel(LabelString);
        ylabel('% of bouts');
        title(BirdNames{i});

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
switch StructString
    case 'FF'
        set(gcf, 'Position', [131 558 400 400]);
    case 'LogAmplitude'
        set(gcf, 'Position', [131 558 1100 400]);
end
p = panel();
switch StructString
    case 'FF'
        p.pack('h', {1});
    case 'LogAmplitude'
        p.pack('h', {1/3 1/3 1/3});
end

p.de.margin = 10;

switch StructString
    case 'FF'
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
        
        
    case 'LogAmplitude'

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
end

set(gca, 'XTick', [1 2], 'XTickLabel', {'L0 Dir' 'Undir'});
axis auto;
ylabel(['Mean of ', LabelString]);

% Stats
switch StructString
    case 'FF'
        Group_L0Dir_UnDir_PValue.Mean = signrank(AllL0Dir_UnDirCVs(1, find(~isnan(AllL0Dir_UnDirCVs(1,:)))), AllL0Dir_UnDirCVs(2, find(~isnan(AllL0Dir_UnDirCVs(1,:)))));
        sigstar({[1 2]}, Group_L0Dir_UnDir_PValue.Mean);
        if (nanmean(AllL0Dir_UnDirCVs(1,:)) > nanmean(AllL0Dir_UnDirCVs(2,:)))
            disp(['Sign rank test comparing group mean CVs for L0 directed and Undirected, p = ', num2str(Group_L0Dir_UnDir_PValue.Mean), '; Dir > Undir']);
        else
            if (nanmean(AllL0Dir_UnDirCVs(1,:)) < nanmean(AllL0Dir_UnDirCVs(2,:)))
                disp(['Sign rank test comparing group mean CVs for L0 directed and Undirected, p = ', num2str(Group_L0Dir_UnDir_PValue.Mean), '; Dir < Undir']);
            else
                disp(['Sign rank test comparing group means for L0 directed and Undirected, p = ', num2str(Group_L0Dir_UnDir_PValue.Mean), '; Dir = Undir']);
            end
        
        end
    case 'LogAmplitude'
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

end

% Now to plot the mean and std across distances for D, D/UN and UN for all
% birds together as percent change relative to undirected.
% This is done only for birds where the means are already different between
% L0 directed and undirected as measured by the earlier p-value

% GroupDataFig = figure;
% set(gcf, 'Color', 'w');
% set(gcf, 'Position', [800 100 500 800]);
% p = panel();
% p.pack({1/3 1/3 1/3});
% p.de.margin = 15;
% p(1).select();
% hold on;
% p(2).select();
% hold on;
% p(3).select();
% hold on;

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
Index = 0;

CV = @(x) std(x)/mean(x);

for i = 1:length(BirdNames),
    for SyllIndex = 1:length(GroupData{i}),
        Index = Index + 1;
        AllDirSongIndices = find(GroupData{i}{SyllIndex}(:,3) == 1);

        % First for directed songs
        for j = 1:6,
            DirSongIndices = find((GroupData{i}{SyllIndex}(:,3) == 1) & (GroupData{i}{SyllIndex}(:,2) == j));
            DirSong_NumTrials(i,j) = length(DirSongIndices);
            if (length(DirSongIndices) >= MinTrialNo)
                DistanceDirMeans(Index, j) = mean(GroupData{i}{SyllIndex}(DirSongIndices,1));
                DistanceDirSEMs(Index, j) = std(GroupData{i}{SyllIndex}(DirSongIndices,1))/sqrt(length(DirSongIndices));
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

        figure(DirUnDirComparisonFig);
        subplot(4, 8, Index);
        hold on;
        switch StructString
            case 'FF'
                plot(Distances, DistanceDirCVs(Index,:), 'rs-', 'LineWidth', 1);
                plot(repmat(Distances(:), 1, 2)', DistanceDirCVCIs{Index}', 'r')

            %     plot(Distances, DistanceDirUnDirMeans(Index,:), 'ks-', 'LineWidth', 1);
            %     plot(repmat(Distances(:), 1, 2)', DistanceDirUnDirCIs{Index}', 'k')
            %     
                plot(Distances, DistanceUnDirCVs(Index,:), 'bs-', 'LineWidth', 1);
                plot(repmat(Distances(:), 1, 2)', DistanceUnDirCVCIs{Index}', 'b');
                plot(repmat(Distances(:),1,2), repmat(DistanceUnDirCVCIs{Index}(end,:), length(Distances), 1), 'k--');

                
            case 'LogAmplitude'
                plot(Distances, DistanceDirMeans(Index,:), 'rs-', 'LineWidth', 1);
                plot(repmat(Distances(:), 1, 2)', DistanceDirCIs{Index}', 'r')

            %     plot(Distances, DistanceDirUnDirMeans(Index,:), 'ks-', 'LineWidth', 1);
            %     plot(repmat(Distances(:), 1, 2)', DistanceDirUnDirCIs{Index}', 'k')
            %     
                plot(Distances, DistanceUnDirMeans(Index,:), 'bs-', 'LineWidth', 1);
                plot(repmat(Distances(:), 1, 2)', DistanceUnDirCIs{Index}', 'b');
                plot(repmat(Distances(:),1,2), repmat(DistanceUnDirCIs{Index}(end,:), length(Distances), 1), 'k--');
        end
        
        xlabel('Distance from female (cm)');
        axis tight;
        Temp = axis;
        Temp = [-0.5 250 0.95*Temp(3) 1.05*Temp(4)];
        axis(Temp);
        ylabel(LabelString);
        title({[BirdNames{i}, ': r=', num2str(Distance_DirSong_Corr{i}(SyllIndex,1)), '; p=', num2str(Distance_DirSong_Corr{i}(SyllIndex,2)), ';p=', num2str(KWP_DirSong{i}(SyllIndex))]; [BirdNames{i}, ': r=', num2str(Distance_UnDirSong_Corr{i}(SyllIndex,1)), '; p=', num2str(Distance_UnDirSong_Corr{i}(SyllIndex,2)), ';p=', num2str(KWP_UnDirSong{i}(SyllIndex))]}, 'Color', 'k', 'FontSize', 8, 'FontWeight', 'normal');
        % legend('D', 'DUN', 'UN');

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
