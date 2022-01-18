function [] = Harini_PlotFFCV(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% == FF analysis for different locations ==================================
% Plots to make
% 1. For each bird
% a. For each syllable, we want distributions in the 3 different conditions
% - all bouts, directed and undirected
% b. We want for each syllable, mean FF and CV across the different
% locations for all 3 conditions
% 2. Group data
% a. Across all syllables
% b. Across all birds
% =========================================================================

Colours = 'rgbcmk';
Symbols = '+o<sd^>*pxvh';
Conditions = [{'L0'} {'L1'} {'L2'} {'L3'} {'L4'} {'UN'}];
MinTrialNo = 5;

for i = 1:length(IndividualBirds),
    NumConditions(i) = length(IndividualBirds(i).ConditionDistances);
end

[MaxNumConditions, MaxNumConditionsIndex] = max(NumConditions);

ConditionDistances = IndividualBirds(MaxNumConditionsIndex).ConditionDistances;

Fid = fopen('FFCV_Correlations.txt', 'w');
fprintf(Fid, 'Bird Name\tPearsons correlation\t\t\t\t\t\tSpearmans correlation\t\t\t\t\t\n');
fprintf(Fid, '\tAllBouts\t\tDirected bouts\t\tUndirected bouts\t\tAllBouts\t\tDirected bouts\t\tUndirected bouts\t\n');
fprintf(Fid, '\tr\tp\tr\tp\tr\tp\tr\tp\tr\tp\tr\tp\n');

% First collect data for all 3 types of bouts listed above
for i = 1:length(IndividualBirds),
    ConditionColumnIndex = find(strcmp('Condition', IndividualBirds(i).BoutStatisticsColumnNames));
    
    IndividualBird_AllBouts{i} = [];
    IndividualBird_DirectedBouts{i} = [];
    IndividualBird_UnDirectedBouts{i} = [];
    
    % First find the syllables that are motif syllables and that are
    % syllables for which we have measured FF
    [CommonSylls, ia, ib] = intersect(IndividualBirds(i).SortedBirdParameters(1).MotifLabels(:), IndividualBirds(i).FF_SyllLabels(:));
    
    FFSyllableIndices = ib;
    
    FFIndex = 1;
    for k = FFSyllableIndices(:)',
        % all the syllables that we have calculated FF for
        IndividualBird_AllBouts{i}{FFIndex} = [];
        IndividualBird_DirectedBouts{i}{FFIndex} = [];
        IndividualBird_UnDirectedBouts{i}{FFIndex} = [];   
        
        for j = 1:1:length(IndividualBirds(i).Conditions),
            % First find all the bouts that correspond to the given condition 
            if (strfind(IndividualBirds(i).Conditions{1}, 'L0'))
                ConditionIndices = find(IndividualBirds(i).BoutStatistics(:,ConditionColumnIndex) == j);
            else
                ConditionIndices = find(IndividualBirds(i).BoutStatistics(:,ConditionColumnIndex) == j+1);
            end
            
            if (length(ConditionIndices) >= MinTrialNo)
                % Add to all bouts
                eval(['AllBouts_', IndividualBirds(i).Conditions{j}, '{', num2str(i), '}(', num2str(FFIndex), ') = ', num2str(nanmean(cell2mat(IndividualBirds(i).FFSyll(k).MeanFF(ConditionIndices)))), ';']);
                eval(['AllBouts_CV_', IndividualBirds(i).Conditions{j}, '{', num2str(i), '}(', num2str(FFIndex), ') = ', num2str(nanstd(cell2mat(IndividualBirds(i).FFSyll(k).MeanFF(ConditionIndices)))/nanmean(cell2mat(IndividualBirds(i).FFSyll(k).MeanFF(ConditionIndices)))), ';']);
                IndividualBird_AllBouts{i}{FFIndex} = [IndividualBird_AllBouts{i}{FFIndex}; [cell2mat(IndividualBirds(i).FFSyll(k).MeanFF(ConditionIndices))' ones(size(cell2mat(IndividualBirds(i).FFSyll(k).MeanFF(ConditionIndices))'))*j]];

                % Find all directed only bouts
                if (isempty(strfind(IndividualBirds(i).Conditions{j}, 'UN')))
                    DirSongIndices = ConditionIndices(find(strcmp('D', IndividualBirds(i).BoutCategorisation(ConditionIndices))));
                else
                    DirSongIndices = ConditionIndices;
                end

                if (length(DirSongIndices) >= MinTrialNo)
                    eval(['DirectedBouts_', IndividualBirds(i).Conditions{j}, '{', num2str(i), '}(', num2str(FFIndex), ') = ', num2str(nanmean(cell2mat(IndividualBirds(i).FFSyll(k).MeanFF(DirSongIndices)))), ';']);
                    eval(['DirectedBouts_CV_', IndividualBirds(i).Conditions{j}, '{', num2str(i), '}(', num2str(FFIndex), ') = ', num2str(nanstd(cell2mat(IndividualBirds(i).FFSyll(k).MeanFF(DirSongIndices)))/nanmean(cell2mat(IndividualBirds(i).FFSyll(k).MeanFF(DirSongIndices)))), ';']);
                    IndividualBird_DirectedBouts{i}{FFIndex} = [IndividualBird_DirectedBouts{i}{FFIndex}; [cell2mat(IndividualBirds(i).FFSyll(k).MeanFF(DirSongIndices))' ones(size(cell2mat(IndividualBirds(i).FFSyll(k).MeanFF(DirSongIndices))'))*j]];
                else
                    eval(['DirectedBouts_', IndividualBirds(i).Conditions{j}, '{', num2str(i), '}(', num2str(FFIndex), ') = ', num2str(NaN), ';']);
                    eval(['DirectedBouts_CV_', IndividualBirds(i).Conditions{j}, '{', num2str(i), '}(', num2str(FFIndex), ') = ', num2str(NaN), ';']);
                end

                % Find all undirected only bouts
                if (isempty(strfind(IndividualBirds(i).Conditions{j}, 'UN')))
                    UnDirSongIndices = ConditionIndices(find(strcmp('UN', IndividualBirds(i).BoutCategorisation(ConditionIndices))));
                else
                    UnDirSongIndices = ConditionIndices;
                end

                if (length(UnDirSongIndices) >= MinTrialNo)
                    eval(['UnDirectedBouts_', IndividualBirds(i).Conditions{j}, '{', num2str(i), '}(', num2str(FFIndex), ') = ', num2str(nanmean(cell2mat(IndividualBirds(i).FFSyll(k).MeanFF(UnDirSongIndices)))), ';']);
                    eval(['UnDirectedBouts_CV_', IndividualBirds(i).Conditions{j}, '{', num2str(i), '}(', num2str(FFIndex), ') = ', num2str(nanstd(cell2mat(IndividualBirds(i).FFSyll(k).MeanFF(UnDirSongIndices)))/nanmean(cell2mat(IndividualBirds(i).FFSyll(k).MeanFF(UnDirSongIndices)))), ';']);
                    IndividualBird_UnDirectedBouts{i}{FFIndex} = [IndividualBird_UnDirectedBouts{i}{FFIndex}; [cell2mat(IndividualBirds(i).FFSyll(k).MeanFF(UnDirSongIndices))' ones(size(cell2mat(IndividualBirds(i).FFSyll(k).MeanFF(UnDirSongIndices))'))*j]];
                else
                    eval(['UnDirectedBouts_', IndividualBirds(i).Conditions{j}, '{', num2str(i), '}(', num2str(FFIndex), ') = ', num2str(NaN), ';']);
                    eval(['UnDirectedBouts_CV_', IndividualBirds(i).Conditions{j}, '{', num2str(i), '}(', num2str(FFIndex), ') = ', num2str(NaN), ';']);
                end
            else
                eval(['AllBouts_', IndividualBirds(i).Conditions{j}, '{', num2str(i), '}(', num2str(FFIndex), ') = ', num2str(NaN), ';']);
                eval(['AllBouts_CV_', IndividualBirds(i).Conditions{j}, '{', num2str(i), '}(', num2str(FFIndex), ') = ', num2str(NaN), ';']);
            end
        end
        
        FFIndex = FFIndex + 1;
    end
    
    for FFSylls = 1:length(IndividualBird_AllBouts{i}),
        % Now to do the stats for All bouts - kruskal-wallis and then
        % multcompare - for means
        [p, tabl, stats] = kruskalwallis(IndividualBird_AllBouts{i}{FFSylls}(:,1), IndividualBird_AllBouts{i}{FFSylls}(:,2), 'off');
        SignificanceMatrix{i}{FFSylls} = ones(length(IndividualBirds(i).Conditions))*NaN;
        if (p < 0.05)
            Comparisons = multcompare(stats, 'display', 'off');
            for k = 1:size(Comparisons,1),
                if (Comparisons(k,end) < 0.05)
                    if ((Comparisons(k,3) < 0) && (Comparisons(k,4) < 0)) 
                        SignificanceMatrix{i}{FFSylls}(Comparisons(k,1), Comparisons(k,2)) = -1;
                    else
                        SignificanceMatrix{i}{FFSylls}(Comparisons(k,1), Comparisons(k,2)) = 1;
                    end
                else
                    SignificanceMatrix{i}{FFSylls}(Comparisons(k,1), Comparisons(k,2)) = 0;
                end
            end
        end
        % Now to do the stats for variance comparison
        [p, stats] = vartestn(IndividualBird_AllBouts{i}{FFSylls}(:,1), IndividualBird_AllBouts{i}{FFSylls}(:,2), 'display', 'off');
        VarSignificance{i}(FFSylls) = p;
    end

%     % ANd then do the correlation linear
%     % First Pearson's
%     Indices = find(IndividualBird_AllBouts{i}(:,2) < (length(IndividualBirds(i).Conditions) - 1));
%     [r, p] = corrcoef(IndividualBird_AllBouts{i}(Indices,1), IndividualBird_AllBouts{i}(Indices,2));
%     Pearsons_R(i) = r(1,2);
%     Pearsons_P(i) = p(1,2);
% 
%     % Also non-parametric Spearmans
%     [r, p] = corr(IndividualBird_AllBouts{i}(Indices,1), ConditionDistances(IndividualBird_AllBouts{i}(Indices,2))', 'type', 'Spearman', 'rows', 'complete');
%     Spearman_R(i) = r;
%     Spearman_P(i) = p;
% 
%     % Now for the directed and undirected bouts, but only for birds
%     % that have enough data from atleast 3 different female positions
% 
%     Directed_Pearsons_R(i) = NaN;
%     Directed_Pearsons_P(i) = NaN;
%     Directed_Spearman_R(i) = NaN;
%     Directed_Spearman_P(i) = NaN;
% 
%     UnDirected_Pearsons_R(i) = NaN;
%     UnDirected_Pearsons_P(i) = NaN;
%     UnDirected_Spearman_R(i) = NaN;
%     UnDirected_Spearman_P(i) = NaN;
% 
%     if (length(unique(IndividualBird_DirectedBouts{i}(:,2))) > 3)
%         Indices = find(IndividualBird_DirectedBouts{i}(:,2) < (length(IndividualBirds(i).Conditions) - 1));
%         [r, p] = corrcoef(IndividualBird_DirectedBouts{i}(Indices,1), IndividualBird_DirectedBouts{i}(Indices,2));
%         Directed_Pearsons_R(i) = r(1,2);
%         Directed_Pearsons_P(i) = p(1,2);
% 
%         % Also non-parametric Spearmans
%         [r, p] = corr(IndividualBird_DirectedBouts{i}(Indices,1), ConditionDistances(IndividualBird_DirectedBouts{i}(Indices,2))', 'type', 'Spearman', 'rows', 'complete');
%         Directed_Spearman_R(i) = r;
%         Directed_Spearman_P(i) = p;
%     end
% 
%     if (length(unique(IndividualBird_UnDirectedBouts{i}(:,2))) > 3)
%         Indices = find(IndividualBird_UnDirectedBouts{i}(:,2) < (length(IndividualBirds(i).Conditions) - 1));
%         [r, p] = corrcoef(IndividualBird_UnDirectedBouts{i}(Indices,1), IndividualBird_UnDirectedBouts{i}(Indices,2));
%         UnDirected_Pearsons_R(i) = r(1,2);
%         UnDirected_Pearsons_P(i) = p(1,2);
% 
%         % Also non-parametric Spearmans
%         [r, p] = corr(IndividualBird_UnDirectedBouts{i}(Indices,1), ConditionDistances(IndividualBird_UnDirectedBouts{i}(Indices,2))', 'type', 'Spearman', 'rows', 'complete');
%         UnDirected_Spearman_R(i) = r;
%         UnDirected_Spearman_P(i) = p;
%     end
%     
%     % Now write all correlation co-efficients and p-values to file
%     fprintf(Fid, '%s\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n', BirdNames{i}, Pearsons_R(i), Pearsons_P(i), Directed_Pearsons_R(i), Directed_Pearsons_P(i), UnDirected_Pearsons_R(i), UnDirected_Pearsons_P(i), Spearman_R(i), Spearman_P(i), Directed_Spearman_R(i), Directed_Spearman_P(i), UnDirected_Spearman_R(i), UnDirected_Spearman_P(i));
    
end
fclose(Fid);

% Put all data together
AllBouts = [];
AllBouts_CV = [];
for i = 1:length(AllBouts_L0),
    
    AllBouts = [AllBouts; [AllBouts_L0{i}(:) AllBouts_L1{i}(:) AllBouts_L2{i}(:) AllBouts_L3{i}(:) AllBouts_L4{i}(:) AllBouts_UN{i}(:)]];
    AllBouts_CV = [AllBouts_CV; [AllBouts_CV_L0{i}(:) AllBouts_CV_L1{i}(:) AllBouts_CV_L2{i}(:) AllBouts_CV_L3{i}(:) AllBouts_CV_L4{i}(:) AllBouts_CV_UN{i}(:)]];
end
AllBouts(find(AllBouts == 0)) = NaN;
AllBouts_CV(find(AllBouts_CV == 0)) = NaN;

DirectedBouts = [];
DirectedBouts_CV = [];
for i = 1:length(DirectedBouts_L0),
    DirectedBouts = [DirectedBouts; [DirectedBouts_L0{i}(:) DirectedBouts_L1{i}(:) DirectedBouts_L2{i}(:) DirectedBouts_L3{i}(:) DirectedBouts_L4{i}(:) DirectedBouts_UN{i}(:)]];
    DirectedBouts_CV = [DirectedBouts_CV; [DirectedBouts_CV_L0{i}(:) DirectedBouts_CV_L1{i}(:) DirectedBouts_CV_L2{i}(:) DirectedBouts_CV_L3{i}(:) DirectedBouts_CV_L4{i}(:) DirectedBouts_CV_UN{i}(:)]];
end
DirectedBouts(find(DirectedBouts == 0)) = NaN;
DirectedBouts_CV(find(DirectedBouts_CV == 0)) = NaN;

UnDirectedBouts = [];
UnDirectedBouts_CV = [];
for i = 1:length(UnDirectedBouts_L0),
    UnDirectedBouts = [UnDirectedBouts; [UnDirectedBouts_L0{i}(:) UnDirectedBouts_L1{i}(:) UnDirectedBouts_L2{i}(:) UnDirectedBouts_L3{i}(:) UnDirectedBouts_L4{i}(:) UnDirectedBouts_UN{i}(:)]];
    UnDirectedBouts_CV = [UnDirectedBouts_CV; [UnDirectedBouts_CV_L0{i}(:) UnDirectedBouts_CV_L1{i}(:) UnDirectedBouts_CV_L2{i}(:) UnDirectedBouts_CV_L3{i}(:) UnDirectedBouts_CV_L4{i}(:) UnDirectedBouts_CV_UN{i}(:)]];
end
UnDirectedBouts(find(UnDirectedBouts == 0)) = NaN;
UnDirectedBouts_CV(find(UnDirectedBouts_CV == 0)) = NaN;

% Now the plots and the statistics

% First the distributions along with the means for all 3 types of bouts in
% 3 subplots of a figure

for i = 1:length(IndividualBirds),
    % Now plot the distributions
    DistributionFig(i) = figure;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [680 166 1000 850]);
    
    MeanCVFig(i) = figure;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [680 166 1000 850]);
    
    for FFSyll = 1:length(IndividualBird_AllBouts{i}),
        figure(DistributionFig(i));
        subplot(length(IndividualBird_AllBouts{i}), 3, (FFSyll - 1)*3 + 1);
        hold on;
        
        PlotLegend{i} = [];
        UniqueConditions = unique(IndividualBird_AllBouts{i}{FFSyll}(:,2));
        
        for j = UniqueConditions(:)',
            ConditionIndices = find(IndividualBird_AllBouts{i}{FFSyll}(:,2) == j);
            plot(min(IndividualBird_AllBouts{i}{FFSyll}(:,1)):5:max(IndividualBird_AllBouts{i}{FFSyll}(:,1)), histc(IndividualBird_AllBouts{i}{FFSyll}(ConditionIndices, 1), min(IndividualBird_AllBouts{i}{FFSyll}(:,1)):5:max(IndividualBird_AllBouts{i}{FFSyll}(:,1)))/length(ConditionIndices), [Colours(mod(j, length(Colours)) + 1), 'o-'], 'LineWidth', 2);
            if (strfind(IndividualBirds(i).Conditions{1}, 'L0'))
                PlotLegend{i}{end+1} = IndividualBirds(i).Conditions{j};
            else
                PlotLegend{i}{end+1} = IndividualBirds(i).Conditions{j-1};
            end
        end

        for j = 1:length(IndividualBirds(i).Conditions),
            ConditionIndices = find(IndividualBird_AllBouts{i}{FFSyll}(:,2) == j);
            plot(mean(IndividualBird_AllBouts{i}{FFSyll}(ConditionIndices,1)), 1 + j*0.05, [Colours(mod(j, length(Colours)) + 1), 's'], 'MarkerSize', 8, 'LineWidth', 1.5);
            plot([(mean(IndividualBird_AllBouts{i}{FFSyll}(ConditionIndices, 1)) - std(IndividualBird_AllBouts{i}{FFSyll}(ConditionIndices, 1))) (mean(IndividualBird_AllBouts{i}{FFSyll}(ConditionIndices, 1)) + std(IndividualBird_AllBouts{i}{FFSyll}(ConditionIndices, 1)))], [(1 + j*0.05) (1 + j*0.05)], Colours(mod(j, length(Colours)) + 1), 'LineWidth', 1.5);
        end

        Temp = [min(IndividualBird_AllBouts{i}{FFSyll}(:,1))-1 max(IndividualBird_AllBouts{i}{FFSyll}(:,1))+15 0 1.35];
        axis(Temp);
        legend(PlotLegend{i});
        xlabel('Syllable FF (Hz)');
        ylabel('Fraction of bouts');
        title([BirdNames{i}, ': All bouts']);

        subplot(length(IndividualBird_AllBouts{i}), 3, (FFSyll - 1)*3 + 2);
        hold on;
        PlotLegend{i} = [];
        UniqueConditions = unique(IndividualBird_DirectedBouts{i}{FFSyll}(:,2));
        
        for j = UniqueConditions(:)',
            ConditionIndices = find(IndividualBird_DirectedBouts{i}{FFSyll}(:,2) == j);
            plot(min(IndividualBird_AllBouts{i}{FFSyll}(:,1)):5:max(IndividualBird_AllBouts{i}{FFSyll}(:,1)), histc(IndividualBird_DirectedBouts{i}{FFSyll}(ConditionIndices, 1), min(IndividualBird_AllBouts{i}{FFSyll}(:,1)):5:max(IndividualBird_AllBouts{i}{FFSyll}(:,1)))/length(ConditionIndices), [Colours(mod(j, length(Colours)) + 1), 'o-'], 'LineWidth', 2);
            if (strfind(IndividualBirds(i).Conditions{1}, 'L0'))
                PlotLegend{i}{end+1} = IndividualBirds(i).Conditions{j};
            else
                PlotLegend{i}{end+1} = IndividualBirds(i).Conditions{j-1};
            end
        end

        for j = 1:length(IndividualBirds(i).Conditions),
            ConditionIndices = find(IndividualBird_DirectedBouts{i}{FFSyll}(:,2) == j);
            plot(mean(IndividualBird_DirectedBouts{i}{FFSyll}(ConditionIndices,1)), 1 + j*0.05, [Colours(mod(j, length(Colours)) + 1), 's'], 'MarkerSize', 8, 'LineWidth', 1.5);
            plot([(mean(IndividualBird_DirectedBouts{i}{FFSyll}(ConditionIndices, 1)) - std(IndividualBird_DirectedBouts{i}{FFSyll}(ConditionIndices, 1))) (mean(IndividualBird_DirectedBouts{i}{FFSyll}(ConditionIndices, 1)) + std(IndividualBird_DirectedBouts{i}{FFSyll}(ConditionIndices, 1)))], [(1 + j*0.05) (1 + j*0.05)], Colours(mod(j, length(Colours)) + 1), 'LineWidth', 1.5);
        end

        Temp = [min(IndividualBird_AllBouts{i}{FFSyll}(:,1))-1 max(IndividualBird_AllBouts{i}{FFSyll}(:,1))+15 0 1.35];
        axis(Temp);
        legend(PlotLegend{i});
        xlabel('Syllable FF (Hz)');
        ylabel('Fraction of bouts');
        title([BirdNames{i}, ': Directed bouts']);

        subplot(length(IndividualBird_AllBouts{i}), 3, (FFSyll - 1)*3 + 3);
        hold on;
        PlotLegend{i} = [];
        UniqueConditions = unique(IndividualBird_UnDirectedBouts{i}{FFSyll}(:,2));
        
        for j = UniqueConditions(:)',
            ConditionIndices = find(IndividualBird_UnDirectedBouts{i}{FFSyll}(:,2) == j);
            plot(min(IndividualBird_AllBouts{i}{FFSyll}(:,1)):5:max(IndividualBird_AllBouts{i}{FFSyll}(:,1)), histc(IndividualBird_UnDirectedBouts{i}{FFSyll}(ConditionIndices, 1), min(IndividualBird_AllBouts{i}{FFSyll}(:,1)):5:max(IndividualBird_AllBouts{i}{FFSyll}(:,1)))/length(ConditionIndices), [Colours(mod(j, length(Colours)) + 1), 'o-'], 'LineWidth', 2);
            if (strfind(IndividualBirds(i).Conditions{1}, 'L0'))
                PlotLegend{i}{end+1} = IndividualBirds(i).Conditions{j};
            else
                PlotLegend{i}{end+1} = IndividualBirds(i).Conditions{j-1};
            end
        end

        for j = 1:length(IndividualBirds(i).Conditions),
            ConditionIndices = find(IndividualBird_UnDirectedBouts{i}{FFSyll}(:,2) == j);
            plot(mean(IndividualBird_UnDirectedBouts{i}{FFSyll}(ConditionIndices,1)), 1 + j*0.05, [Colours(mod(j, length(Colours)) + 1), 's'], 'MarkerSize', 8, 'LineWidth', 1.5);
            plot([(mean(IndividualBird_UnDirectedBouts{i}{FFSyll}(ConditionIndices, 1)) - std(IndividualBird_UnDirectedBouts{i}{FFSyll}(ConditionIndices, 1))) (mean(IndividualBird_UnDirectedBouts{i}{FFSyll}(ConditionIndices, 1)) + std(IndividualBird_UnDirectedBouts{i}{FFSyll}(ConditionIndices, 1)))], [(1 + j*0.05) (1 + j*0.05)], Colours(mod(j, length(Colours)) + 1), 'LineWidth', 1.5);
        end

        Temp = [min(IndividualBird_AllBouts{i}{FFSyll}(:,1))-1 max(IndividualBird_AllBouts{i}{FFSyll}(:,1))+15 0 1.35];
        axis(Temp);
        legend(PlotLegend{i});
        xlabel('Syllable FF (Hz)');
        ylabel('Fraction of bouts');
        title([BirdNames{i}, ': Undirected bouts']);

        figure(MeanCVFig(i));
        subplot(length(IndividualBird_AllBouts{i}), 3, (FFSyll - 1)*3 + 1);
        hold on;
        MeanSEMValues = [];
        CVValues = [];
        for j = min(IndividualBird_AllBouts{i}{FFSyll}(:,2)):1:max(IndividualBird_AllBouts{i}{FFSyll}(:,2)),
            Indices = find(IndividualBird_AllBouts{i}{FFSyll}(:,2) == j);
            if (~isempty(Indices))
                MeanSEMValues(end+1,:) = [j nanmean(IndividualBird_AllBouts{i}{FFSyll}(Indices,1)) nanstd(IndividualBird_AllBouts{i}{FFSyll}(Indices,1))/sqrt(length(Indices))];
                CVValues(end+1,:) = [j nanstd(IndividualBird_AllBouts{i}{FFSyll}(Indices,1))/nanmean(IndividualBird_AllBouts{i}{FFSyll}(Indices,1))];
            else
                MeanSEMValues(end+1,:) = [j NaN NaN];
                CVValues(end+1,:) = [j NaN];
            end
        end
        errorbar(ConditionDistances(MeanSEMValues(:,1)), MeanSEMValues(:,2), MeanSEMValues(:,3), 'ko-', 'LineWidth', 2, 'MarkerSize', 7);
        
        subplot(length(IndividualBird_AllBouts{i}), 3, (FFSyll - 1)*3 + 2);
        hold on;
        plot(ConditionDistances(CVValues(:,1)), CVValues(:,2), 'ko-', 'LineWidth', 2, 'MarkerSize', 7);
        
        MeanSEMValues = [];
        CVValues = [];
        for j = min(IndividualBird_DirectedBouts{i}{FFSyll}(:,2)):1:max(IndividualBird_DirectedBouts{i}{FFSyll}(:,2)),
            Indices = find(IndividualBird_DirectedBouts{i}{FFSyll}(:,2) == j);
            if (~isempty(Indices))
                MeanSEMValues(end+1,:) = [j nanmean(IndividualBird_DirectedBouts{i}{FFSyll}(Indices,1)) nanstd(IndividualBird_DirectedBouts{i}{FFSyll}(Indices,1))/sqrt(length(Indices))];
                CVValues(end+1,:) = [j nanstd(IndividualBird_DirectedBouts{i}{FFSyll}(Indices,1))/nanmean(IndividualBird_DirectedBouts{i}{FFSyll}(Indices,1))];
            else
                MeanSEMValues(end+1,:) = [j NaN NaN];
                CVValues(end+1,:) = [j NaN];
            end
        end
        subplot(length(IndividualBird_AllBouts{i}), 3, (FFSyll - 1)*3 + 1);
        errorbar(ConditionDistances(MeanSEMValues(:,1)), MeanSEMValues(:,2), MeanSEMValues(:,3), 'ro-', 'LineWidth', 2, 'MarkerSize', 7);
        subplot(length(IndividualBird_AllBouts{i}), 3, (FFSyll - 1)*3 + 2);
        plot(ConditionDistances(CVValues(:,1)), CVValues(:,2), 'ro-', 'LineWidth', 2, 'MarkerSize', 7);
        
        MeanSEMValues = [];
        CVValues = [];
        for j = min(IndividualBird_UnDirectedBouts{i}{FFSyll}(:,2)):1:max(IndividualBird_UnDirectedBouts{i}{FFSyll}(:,2)),
            Indices = find(IndividualBird_UnDirectedBouts{i}{FFSyll}(:,2) == j);
            if (~isempty(Indices))
                MeanSEMValues(end+1,:) = [j nanmean(IndividualBird_UnDirectedBouts{i}{FFSyll}(Indices,1)) nanstd(IndividualBird_UnDirectedBouts{i}{FFSyll}(Indices,1))/sqrt(length(Indices))];
                CVValues(end+1,:) = [j nanstd(IndividualBird_UnDirectedBouts{i}{FFSyll}(Indices,1))/nanmean(IndividualBird_UnDirectedBouts{i}{FFSyll}(Indices,1))];
            else
                MeanSEMValues(end+1,:) = [j NaN NaN];
                CVValues(end+1,:) = [j NaN];
            end
        end
        subplot(length(IndividualBird_AllBouts{i}), 3, (FFSyll - 1)*3 + 1);
        errorbar(ConditionDistances(MeanSEMValues(:,1)), MeanSEMValues(:,2), MeanSEMValues(:,3), 'bo-', 'LineWidth', 2, 'MarkerSize', 7);
        subplot(length(IndividualBird_AllBouts{i}), 3, (FFSyll - 1)*3 + 2);
        plot(ConditionDistances(CVValues(:,1)), CVValues(:,2), 'bo-', 'LineWidth', 2, 'MarkerSize', 7);

        axis tight;
        Temp = axis;
        Temp = [(ConditionDistances(1) - 1) (ConditionDistances(end) + 1) 0.98*Temp(3) 1.02*Temp(4)];
        axis(Temp);
        legend('All bouts', 'Directed bouts', 'Undirected bouts', 'Location', 'southeast');
        ylabel('Mean Syllable FF (Hz)');
        xlabel('Distance of female (cm)');
        title([BirdNames{i}]);

        subplot(length(IndividualBird_AllBouts{i}), 3, (FFSyll - 1)*3 + 2);
        axis tight;
        Temp = axis;
        Temp = [(ConditionDistances(1) - 1) (ConditionDistances(end) + 1) 0.98*Temp(3) 1.02*Temp(4)];
        axis(Temp);
        legend('All bouts', 'Directed bouts', 'Undirected bouts', 'Location', 'southeast');
        ylabel('CV of Syllable FF');
        xlabel('Distance of female (cm)');
        title([BirdNames{i}]);
        
        subplot(length(IndividualBird_AllBouts{i}), 3, (FFSyll - 1)*3 + 3);
        TempImage = imagesc(SignificanceMatrix{i}{FFSyll});
        set(TempImage, 'AlphaData', ~isnan(SignificanceMatrix{i}{FFSyll}));
        if (min(SignificanceMatrix{i}{FFSyll}(:)) < 0)
            colormap([0 0 1; 0.5 0.5 0.5; 1 0 0]);
        else
            colormap([0.5 0.5 0.5; 1 0 0]);
        end
        colorbar;
        for j = 1:length(IndividualBirds(i).Conditions),
            if (isempty(strfind(IndividualBirds(i).Conditions{j}, 'UN')))
                AxisLabelString{j} = [num2str(IndividualBirds(i).ConditionDistances(j)), ' cm'];
            else
                AxisLabelString{j} = IndividualBirds(i).Conditions{j};
            end
        end
        set(gca, 'XTick', [1:1:length(IndividualBirds(i).Conditions)], 'XTickLabel', AxisLabelString);
        set(gca, 'YTick', [1:1:length(IndividualBirds(i).Conditions)], 'YTickLabel', AxisLabelString);
    end
end


% Normalise as % of undirected and also as delta change relative to
% undirected
for i = 1:size(AllBouts,1),
    PercentChange_AllBouts(i,:) = ((AllBouts(i,:) - AllBouts(i,end)) * 100)/AllBouts(i,end);
    DeltaChange_AllBouts(i,:) = (AllBouts(i,:) - AllBouts(i,end));
end

% First for percent change
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 289 750 675]);
for i = 1:size(AllBouts,2),
    AllBout_Bar(i) = bar(ConditionDistances(i), nanmean(PercentChange_AllBouts(:,i)));
    set(AllBout_Bar(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'BarWidth', 10, 'LineWidth', 2);
    hold on;
    errorbar(ConditionDistances(i), nanmean(PercentChange_AllBouts(:,i)), nanstd(PercentChange_AllBouts(:,i))/sqrt(length(find(~isnan(PercentChange_AllBouts(:,i))))), 'ks-', 'MarkerSize', 10, 'LineWidth', 2);
end
plot(ConditionDistances, nanmean(PercentChange_AllBouts), 'k', 'LineWidth', 2);

for i = 1:size(AllBouts,1),
    PercentChangeLine(i) = plot(ConditionDistances, PercentChange_AllBouts(i,:), [Symbols(mod(i,length(Symbols))+1), Colours(ceil(i/length(Symbols))), '-'], 'LineWidth', 0.5, 'MarkerSize', 5);
end
legend(PercentChangeLine, BirdNames);

clear AxisLabelString;
for j = 1:length(ConditionDistances),
    if (isempty(strfind(IndividualBirds(MaxNumConditionsIndex).Conditions{j}, 'UN')))
        AxisLabelString{j} = num2str(ConditionDistances(j));
    else
        AxisLabelString{j} = 'UN';
    end
end

set(gca, 'XTick', ConditionDistances, 'XTickLabel', AxisLabelString);
xlabel('Distance of female cage (cm)');
ylabel('% change in mean FF relative to Undirected song');
title(['All bouts: (n=', num2str(length(BirdNames)), ') birds']);
axis tight;
Temp = axis;
Temp = [(ConditionDistances(1) - 12) (ConditionDistances(end) + 12) 1.02*Temp(3) 1.2*Temp(4)];
axis(Temp);

ConditionIndexMatrix = repmat(1:1:size(PercentChange_AllBouts,2), size(PercentChange_AllBouts, 1), 1);
% Omit nan rows for doing kruskalwallis
[Nanrows, Nancols] = find(isnan(PercentChange_AllBouts));
RowsToUse = setdiff(1:1:size(PercentChange_AllBouts,1), unique(Nanrows));

DataMatrix = PercentChange_AllBouts(RowsToUse,:);
IndexMatrix = ConditionIndexMatrix(RowsToUse,:);
[p, tabl, stats] = kruskalwallis(DataMatrix(:), IndexMatrix(:), 'off');
if (p < 0.05)
    Comparisons = multcompare(stats, 'display', 'off');
    % only plotting * for anything that is significantly different from
    % undirected
    SigIndices = find((Comparisons(:,end) < 0.05) & (Comparisons(:,2) == length(ConditionDistances)));
    for j = 1:length(SigIndices),
        plot(ConditionDistances(Comparisons(SigIndices(j,1))), 1.15*(Temp(4))/1.2, 'k*', 'MarkerSize', 10, 'LineWidth', 2);
    end
end

% Now calculate correlation to the directed data alone and for the birds that have
% all distances
DataMatrix = DataMatrix(1:(length(ConditionDistances)-1),:);
IndexMatrix = IndexMatrix(1:(length(ConditionDistances) - 1), :);
[AllBird_PearsonR, AllBird_PearsonP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:));
[AllBird_SpearmanR, AllBirds_SpearmanP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:), 'type', 'Spearman', 'rows', 'complete');
text(30, Temp(4)/1.2, ['r = ', num2str(AllBird_PearsonR), '; p = ', num2str(AllBird_PearsonP)], 'FontSize', 12, 'FontWeight', 'bold');

% First for delta change
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 289 750 675]);
for i = 1:size(AllBouts,2),
    AllBout_Bar(i) = bar(ConditionDistances(i), nanmean(DeltaChange_AllBouts(:,i)));
    set(AllBout_Bar(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'BarWidth', 10, 'LineWidth', 2);
    hold on;
    errorbar(ConditionDistances(i), nanmean(DeltaChange_AllBouts(:,i)), nanstd(DeltaChange_AllBouts(:,i))/sqrt(length(find(~isnan(DeltaChange_AllBouts(:,i))))), 'ks-', 'MarkerSize', 10, 'LineWidth', 2);
end
plot(ConditionDistances, nanmean(DeltaChange_AllBouts), 'k', 'LineWidth', 2);

for i = 1:size(AllBouts,1),
    DeltaChangeLine(i) = plot(ConditionDistances, DeltaChange_AllBouts(i,:), [Symbols(mod(i,length(Symbols))+1), Colours(ceil(i/length(Symbols))), '-'], 'LineWidth', 0.5, 'MarkerSize', 5);
end
legend(DeltaChangeLine, BirdNames);

clear AxisLabelString;
for j = 1:length(ConditionDistances),
    if (isempty(strfind(IndividualBirds(MaxNumConditionsIndex).Conditions{j}, 'UN')))
        AxisLabelString{j} = num2str(ConditionDistances(j));
    else
        AxisLabelString{j} = 'UN';
    end
end

set(gca, 'XTick', ConditionDistances, 'XTickLabel', AxisLabelString);
xlabel('Distance of female cage (cm)');
ylabel('\Delta change in mean FF relative to Undirected song');
title(['All bouts: (n=', num2str(length(BirdNames)), ') birds']);
axis tight;
Temp = axis;
Temp = [(ConditionDistances(1) - 12) (ConditionDistances(end) + 12) 1.02*Temp(3) 1.2*Temp(4)];
axis(Temp);

ConditionIndexMatrix = repmat(1:1:size(DeltaChange_AllBouts,2), size(DeltaChange_AllBouts, 1), 1);
% Omit nan rows for doing kruskalwallis
[Nanrows, Nancols] = find(isnan(DeltaChange_AllBouts));
RowsToUse = setdiff(1:1:size(DeltaChange_AllBouts,1), unique(Nanrows));

DataMatrix = DeltaChange_AllBouts(RowsToUse,:);
IndexMatrix = ConditionIndexMatrix(RowsToUse,:);
[p, tabl, stats] = kruskalwallis(DataMatrix(:), IndexMatrix(:), 'off');
if (p < 0.05)
    Comparisons = multcompare(stats, 'display', 'off');
    % only plotting * for anything that is significantly different from
    % undirected
    SigIndices = find((Comparisons(:,end) < 0.05) & (Comparisons(:,2) == length(ConditionDistances)));
    for j = 1:length(SigIndices),
        plot(ConditionDistances(Comparisons(SigIndices(j,1))), 1.15*(Temp(4))/1.2, 'k*', 'MarkerSize', 10, 'LineWidth', 2);
    end
end

% Now calculate correlation to the directed data alone and for the birds that have
% all distances
DataMatrix = DataMatrix(1:(length(ConditionDistances)-1),:);
IndexMatrix = IndexMatrix(1:(length(ConditionDistances) - 1), :);
[AllBird_PearsonR, AllBird_PearsonP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:));
[AllBird_SpearmanR, AllBirds_SpearmanP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:), 'type', 'Spearman', 'rows', 'complete');
text(30, Temp(4)/1.2, ['r = ', num2str(AllBird_PearsonR), '; p = ', num2str(AllBird_PearsonP)], 'FontSize', 12, 'FontWeight', 'bold');

% Now to do the same thing for all directed only bouts

% Normalise as % of undirected and also as delta change relative to
% undirected
for i = 1:size(DirectedBouts,1),
    PercentChange_DirectedBouts(i,:) = ((DirectedBouts(i,:) - DirectedBouts(i,end)) * 100)/DirectedBouts(i,end);
    DeltaChange_DirectedBouts(i,:) = (DirectedBouts(i,:) - DirectedBouts(i,end));
end

% First for percent change
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 289 750 675]);
for i = 1:size(DirectedBouts,2),
    AllBout_Bar(i) = bar(ConditionDistances(i), nanmean(PercentChange_DirectedBouts(:,i)));
    set(AllBout_Bar(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'BarWidth', 10, 'LineWidth', 2);
    hold on;
    errorbar(ConditionDistances(i), nanmean(PercentChange_DirectedBouts(:,i)), nanstd(PercentChange_DirectedBouts(:,i))/sqrt(length(find(~isnan(PercentChange_DirectedBouts(:,i))))), 'ks-', 'MarkerSize', 10, 'LineWidth', 2);
end
plot(ConditionDistances, nanmean(PercentChange_DirectedBouts), 'k', 'LineWidth', 2);

for i = 1:size(DirectedBouts,1),
    PercentChangeLine(i) = plot(ConditionDistances, PercentChange_DirectedBouts(i,:), [Symbols(mod(i,length(Symbols))+1), Colours(ceil(i/length(Symbols))), '-'], 'LineWidth', 0.5, 'MarkerSize', 5);
end
legend(PercentChangeLine, BirdNames);

clear AxisLabelString;
for j = 1:length(ConditionDistances),
    if (isempty(strfind(IndividualBirds(MaxNumConditionsIndex).Conditions{j}, 'UN')))
        AxisLabelString{j} = num2str(ConditionDistances(j));
    else
        AxisLabelString{j} = 'UN';
    end
end

set(gca, 'XTick', ConditionDistances, 'XTickLabel', AxisLabelString);
xlabel('Distance of female cage (cm)');
ylabel('% change in mean FF relative to Undirected song');
title(['Directed bouts: (n=', num2str(length(BirdNames)), ') birds']);
axis tight;
Temp = axis;
Temp = [(ConditionDistances(1) - 12) (ConditionDistances(end) + 12) 1.02*Temp(3) 1.2*Temp(4)];
axis(Temp);

ConditionIndexMatrix = repmat(1:1:size(PercentChange_DirectedBouts,2), size(PercentChange_DirectedBouts, 1), 1);
% Omit nans for doing kruskalwallis
%DataMatrix = PercentChange_DirectedBouts(RowsToUse,:);
IndexMatrix = ConditionIndexMatrix;
NotNaNIndices = find(~isnan(DataMatrix));
[p, tabl, stats] = kruskalwallis(DataMatrix(NotNaNIndices), IndexMatrix(NotNaNIndices), 'off');
if (p < 0.05)
    Comparisons = multcompare(stats, 'display', 'off');
    % only plotting * for anything that is significantly different from
    % undirected
    SigIndices = find((Comparisons(:,end) < 0.05) & (Comparisons(:,2) == length(ConditionDistances)));
    for j = 1:length(SigIndices),
        plot(ConditionDistances(Comparisons(SigIndices(j,1))), 1.15*(Temp(4))/1.2, 'k*', 'MarkerSize', 10, 'LineWidth', 2);
    end
end

% Now calculate correlation to the directed data alone and for the birds that have
% all distances
% Take only rows that have at least 3 not-nan values in the first 3 columns

[Nanrows, Nancols] = find(isnan(PercentChange_DirectedBouts(:,1:3)));
RowsToUse = setdiff(1:1:size(PercentChange_DirectedBouts,1), unique(Nanrows));

DataMatrix = DeltaChange_AllBouts(RowsToUse,1:3);
IndexMatrix = ConditionIndexMatrix(RowsToUse,1:3);

[AllBird_DirectedPearsonR, AllBird_DirectedPearsonP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:));
[AllBird_DirectedSpearmanR, AllBirds_DirectedSpearmanP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:), 'type', 'Spearman', 'rows', 'complete');
text(30, Temp(4)/1.2, ['r = ', num2str(AllBird_DirectedPearsonR), '; p = ', num2str(AllBird_DirectedPearsonP)], 'FontSize', 12, 'FontWeight', 'bold');


% Now to do the same thing for all undirected only bouts

% Normalise as % of undirected and also as delta change relative to
% undirected
for i = 1:size(UnDirectedBouts,1),
    PercentChange_UnDirectedBouts(i,:) = ((UnDirectedBouts(i,:) - UnDirectedBouts(i,end)) * 100)/UnDirectedBouts(i,end);
    DeltaChange_UnDirectedBouts(i,:) = (UnDirectedBouts(i,:) - UnDirectedBouts(i,end));
end

% First for percent change
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 289 750 675]);
for i = 1:size(UnDirectedBouts,2),
    AllBout_Bar(i) = bar(ConditionDistances(i), nanmean(PercentChange_UnDirectedBouts(:,i)));
    set(AllBout_Bar(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'BarWidth', 10, 'LineWidth', 2);
    hold on;
    errorbar(ConditionDistances(i), nanmean(PercentChange_UnDirectedBouts(:,i)), nanstd(PercentChange_UnDirectedBouts(:,i))/sqrt(length(find(~isnan(PercentChange_UnDirectedBouts(:,i))))), 'ks-', 'MarkerSize', 10, 'LineWidth', 2);
end
plot(ConditionDistances, nanmean(PercentChange_UnDirectedBouts), 'k', 'LineWidth', 2);

for i = 1:size(UnDirectedBouts,1),
    PercentChangeLine(i) = plot(ConditionDistances, PercentChange_UnDirectedBouts(i,:), [Symbols(mod(i,length(Symbols))+1), Colours(ceil(i/length(Symbols))), '-'], 'LineWidth', 0.5, 'MarkerSize', 5);
end
legend(PercentChangeLine, BirdNames);

clear AxisLabelString;
for j = 1:length(ConditionDistances),
    if (isempty(strfind(IndividualBirds(MaxNumConditionsIndex).Conditions{j}, 'UN')))
        AxisLabelString{j} = num2str(ConditionDistances(j));
    else
        AxisLabelString{j} = 'UN';
    end
end

set(gca, 'XTick', ConditionDistances, 'XTickLabel', AxisLabelString);
xlabel('Distance of female cage (cm)');
ylabel('% change in mean FF relative to Undirected song');
title(['UnDirected bouts: (n=', num2str(length(BirdNames)), ') birds']);
axis tight;
Temp = axis;
Temp = [(ConditionDistances(1) - 12) (ConditionDistances(end) + 12) 1.02*Temp(3) 1.2*Temp(4)];
axis(Temp);

ConditionIndexMatrix = repmat(1:1:size(PercentChange_UnDirectedBouts,2), size(PercentChange_UnDirectedBouts, 1), 1);
% Omit nans for doing kruskalwallis
%DataMatrix = PercentChange_UnDirectedBouts(RowsToUse,:);
IndexMatrix = ConditionIndexMatrix;
NotNaNIndices = find(~isnan(DataMatrix));
[p, tabl, stats] = kruskalwallis(DataMatrix(NotNaNIndices), IndexMatrix(NotNaNIndices), 'off');
if (p < 0.05)
    Comparisons = multcompare(stats, 'display', 'off');
    % only plotting * for anything that is significantly different from
    % undirected
    SigIndices = find((Comparisons(:,end) < 0.05) & (Comparisons(:,2) == length(ConditionDistances)));
    for j = 1:length(SigIndices),
        plot(ConditionDistances(Comparisons(SigIndices(j,1))), 1.15*(Temp(4))/1.2, 'k*', 'MarkerSize', 10, 'LineWidth', 2);
    end
end

% Now calculate correlation to the directed data alone and for the birds that have
% all distances
% Take only rows that have at least 3 not-nan values in the last 3 female
% presentation columns

[Nanrows, Nancols] = find(isnan(PercentChange_UnDirectedBouts(:,(length(ConditionDistances) - 3):(length(ConditionDistances) - 1))));
RowsToUse = setdiff(1:1:size(PercentChange_UnDirectedBouts,1), unique(Nanrows));

DataMatrix = DeltaChange_AllBouts(RowsToUse,(length(ConditionDistances) - 3):(length(ConditionDistances) - 1));
IndexMatrix = ConditionIndexMatrix(RowsToUse,(length(ConditionDistances) - 3):(length(ConditionDistances) - 1));

[AllBird_DirectedPearsonR, AllBird_DirectedPearsonP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:));
[AllBird_DirectedSpearmanR, AllBirds_DirectedSpearmanP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:), 'type', 'Spearman', 'rows', 'complete');
text(30, Temp(4)/1.2, ['r = ', num2str(AllBird_DirectedPearsonR), '; p = ', num2str(AllBird_DirectedPearsonP)], 'FontSize', 12, 'FontWeight', 'bold');


% Now do the same for CV
% Normalise as % of undirected and also as delta change relative to
% undirected
for i = 1:size(AllBouts_CV,1),
    PercentChange_AllBouts_CV(i,:) = ((AllBouts_CV(i,:) - AllBouts_CV(i,end)) * 100)/AllBouts_CV(i,end);
    DeltaChange_AllBouts_CV(i,:) = (AllBouts_CV(i,:) - AllBouts_CV(i,end));
end

% First for percent change
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 289 750 675]);
for i = 1:size(AllBouts_CV,2),
    AllBout_Bar(i) = bar(ConditionDistances(i), nanmean(PercentChange_AllBouts_CV(:,i)));
    set(AllBout_Bar(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'BarWidth', 10, 'LineWidth', 2);
    hold on;
    errorbar(ConditionDistances(i), nanmean(PercentChange_AllBouts_CV(:,i)), nanstd(PercentChange_AllBouts_CV(:,i))/sqrt(length(find(~isnan(PercentChange_AllBouts_CV(:,i))))), 'ks-', 'MarkerSize', 10, 'LineWidth', 2);
end
plot(ConditionDistances, nanmean(PercentChange_AllBouts_CV), 'k', 'LineWidth', 2);

for i = 1:size(AllBouts_CV,1),
    PercentChangeLine(i) = plot(ConditionDistances, PercentChange_AllBouts_CV(i,:), [Symbols(mod(i,length(Symbols))+1), Colours(ceil(i/length(Symbols))), '-'], 'LineWidth', 0.5, 'MarkerSize', 5);
end
legend(PercentChangeLine, BirdNames);

clear AxisLabelString;
for j = 1:length(ConditionDistances),
    if (isempty(strfind(IndividualBirds(MaxNumConditionsIndex).Conditions{j}, 'UN')))
        AxisLabelString{j} = num2str(ConditionDistances(j));
    else
        AxisLabelString{j} = 'UN';
    end
end

set(gca, 'XTick', ConditionDistances, 'XTickLabel', AxisLabelString);
xlabel('Distance of female cage (cm)');
ylabel('% change in CV of FF relative to Undirected song');
title(['All bouts: (n=', num2str(length(BirdNames)), ') birds']);
axis tight;
Temp = axis;
Temp = [(ConditionDistances(1) - 12) (ConditionDistances(end) + 12) 1.02*Temp(3) 1.2*Temp(4)];
axis(Temp);

ConditionIndexMatrix = repmat(1:1:size(PercentChange_AllBouts_CV,2), size(PercentChange_AllBouts_CV, 1), 1);
% Omit nan rows for doing kruskalwallis
[Nanrows, Nancols] = find(isnan(PercentChange_AllBouts_CV));
RowsToUse = setdiff(1:1:size(PercentChange_AllBouts_CV,1), unique(Nanrows));

DataMatrix = PercentChange_AllBouts_CV(RowsToUse,:);
IndexMatrix = ConditionIndexMatrix(RowsToUse,:);
[p, tabl, stats] = kruskalwallis(DataMatrix(:), IndexMatrix(:), 'off');
if (p < 0.05)
    Comparisons = multcompare(stats, 'display', 'off');
    % only plotting * for anything that is significantly different from
    % undirected
    SigIndices = find((Comparisons(:,end) < 0.05) & (Comparisons(:,2) == length(ConditionDistances)));
    for j = 1:length(SigIndices),
        plot(ConditionDistances(Comparisons(SigIndices(j,1))), 1.15*(Temp(4))/1.2, 'k*', 'MarkerSize', 10, 'LineWidth', 2);
    end
end

% Now calculate correlation to the directed data alone and for the birds that have
% all distances
DataMatrix = DataMatrix(1:(length(ConditionDistances)-1),:);
IndexMatrix = IndexMatrix(1:(length(ConditionDistances) - 1), :);
[AllBird_PearsonR, AllBird_PearsonP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:));
[AllBird_SpearmanR, AllBirds_SpearmanP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:), 'type', 'Spearman', 'rows', 'complete');
text(30, Temp(4)/1.2, ['r = ', num2str(AllBird_PearsonR), '; p = ', num2str(AllBird_PearsonP)], 'FontSize', 12, 'FontWeight', 'bold');

% First for delta change
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 289 750 675]);
for i = 1:size(AllBouts_CV,2),
    AllBout_Bar(i) = bar(ConditionDistances(i), nanmean(DeltaChange_AllBouts_CV(:,i)));
    set(AllBout_Bar(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'BarWidth', 10, 'LineWidth', 2);
    hold on;
    errorbar(ConditionDistances(i), nanmean(DeltaChange_AllBouts_CV(:,i)), nanstd(DeltaChange_AllBouts_CV(:,i))/sqrt(length(find(~isnan(DeltaChange_AllBouts_CV(:,i))))), 'ks-', 'MarkerSize', 10, 'LineWidth', 2);
end
plot(ConditionDistances, nanmean(DeltaChange_AllBouts_CV), 'k', 'LineWidth', 2);

for i = 1:size(AllBouts_CV,1),
    DeltaChangeLine(i) = plot(ConditionDistances, DeltaChange_AllBouts_CV(i,:), [Symbols(mod(i,length(Symbols))+1), Colours(ceil(i/length(Symbols))), '-'], 'LineWidth', 0.5, 'MarkerSize', 5);
end
legend(DeltaChangeLine, BirdNames);

clear AxisLabelString;
for j = 1:length(ConditionDistances),
    if (isempty(strfind(IndividualBirds(MaxNumConditionsIndex).Conditions{j}, 'UN')))
        AxisLabelString{j} = num2str(ConditionDistances(j));
    else
        AxisLabelString{j} = 'UN';
    end
end

set(gca, 'XTick', ConditionDistances, 'XTickLabel', AxisLabelString);
xlabel('Distance of female cage (cm)');
ylabel('\Delta change in CV of FF relative to Undirected song');
title(['All bouts: (n=', num2str(length(BirdNames)), ') birds']);
axis tight;
Temp = axis;
Temp = [(ConditionDistances(1) - 12) (ConditionDistances(end) + 12) 1.02*Temp(3) 1.2*Temp(4)];
axis(Temp);

ConditionIndexMatrix = repmat(1:1:size(DeltaChange_AllBouts_CV,2), size(DeltaChange_AllBouts_CV, 1), 1);
% Omit nan rows for doing kruskalwallis
[Nanrows, Nancols] = find(isnan(DeltaChange_AllBouts_CV));
RowsToUse = setdiff(1:1:size(DeltaChange_AllBouts_CV,1), unique(Nanrows));

DataMatrix = DeltaChange_AllBouts_CV(RowsToUse,:);
IndexMatrix = ConditionIndexMatrix(RowsToUse,:);
[p, tabl, stats] = kruskalwallis(DataMatrix(:), IndexMatrix(:), 'off');
if (p < 0.05)
    Comparisons = multcompare(stats, 'display', 'off');
    % only plotting * for anything that is significantly different from
    % undirected
    SigIndices = find((Comparisons(:,end) < 0.05) & (Comparisons(:,2) == length(ConditionDistances)));
    for j = 1:length(SigIndices),
        plot(ConditionDistances(Comparisons(SigIndices(j,1))), 1.15*(Temp(4))/1.2, 'k*', 'MarkerSize', 10, 'LineWidth', 2);
    end
end

% Now calculate correlation to the directed data alone and for the birds that have
% all distances
DataMatrix = DataMatrix(1:(length(ConditionDistances)-1),:);
IndexMatrix = IndexMatrix(1:(length(ConditionDistances) - 1), :);
[AllBird_PearsonR, AllBird_PearsonP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:));
[AllBird_SpearmanR, AllBirds_SpearmanP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:), 'type', 'Spearman', 'rows', 'complete');
text(30, Temp(4)/1.2, ['r = ', num2str(AllBird_PearsonR), '; p = ', num2str(AllBird_PearsonP)], 'FontSize', 12, 'FontWeight', 'bold');

% Now to do the same thing for all directed only bouts

% Normalise as % of undirected and also as delta change relative to
% undirected
for i = 1:size(DirectedBouts_CV,1),
    PercentChange_DirectedBouts_CV(i,:) = ((DirectedBouts_CV(i,:) - DirectedBouts_CV(i,end)) * 100)/DirectedBouts_CV(i,end);
    DeltaChange_DirectedBouts_CV(i,:) = (DirectedBouts_CV(i,:) - DirectedBouts_CV(i,end));
end

% First for percent change
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 289 750 675]);
for i = 1:size(DirectedBouts_CV,2),
    AllBout_Bar(i) = bar(ConditionDistances(i), nanmean(PercentChange_DirectedBouts_CV(:,i)));
    set(AllBout_Bar(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'BarWidth', 10, 'LineWidth', 2);
    hold on;
    errorbar(ConditionDistances(i), nanmean(PercentChange_DirectedBouts_CV(:,i)), nanstd(PercentChange_DirectedBouts_CV(:,i))/sqrt(length(find(~isnan(PercentChange_DirectedBouts_CV(:,i))))), 'ks-', 'MarkerSize', 10, 'LineWidth', 2);
end
plot(ConditionDistances, nanmean(PercentChange_DirectedBouts_CV), 'k', 'LineWidth', 2);

for i = 1:size(DirectedBouts_CV,1),
    PercentChangeLine(i) = plot(ConditionDistances, PercentChange_DirectedBouts_CV(i,:), [Symbols(mod(i,length(Symbols))+1), Colours(ceil(i/length(Symbols))), '-'], 'LineWidth', 0.5, 'MarkerSize', 5);
end
legend(PercentChangeLine, BirdNames);

clear AxisLabelString;
for j = 1:length(ConditionDistances),
    if (isempty(strfind(IndividualBirds(MaxNumConditionsIndex).Conditions{j}, 'UN')))
        AxisLabelString{j} = num2str(ConditionDistances(j));
    else
        AxisLabelString{j} = 'UN';
    end
end

set(gca, 'XTick', ConditionDistances, 'XTickLabel', AxisLabelString);
xlabel('Distance of female cage (cm)');
ylabel('% change in CV of FF relative to Undirected song');
title(['Directed bouts: (n=', num2str(length(BirdNames)), ') birds']);
axis tight;
Temp = axis;
Temp = [(ConditionDistances(1) - 12) (ConditionDistances(end) + 12) 1.02*Temp(3) 1.2*Temp(4)];
axis(Temp);

ConditionIndexMatrix = repmat(1:1:size(PercentChange_DirectedBouts_CV,2), size(PercentChange_DirectedBouts_CV, 1), 1);
% Omit nans for doing kruskalwallis
%DataMatrix = PercentChange_DirectedBouts_CV(RowsToUse,:);
IndexMatrix = ConditionIndexMatrix;
NotNaNIndices = find(~isnan(DataMatrix));
[p, tabl, stats] = kruskalwallis(DataMatrix(NotNaNIndices), IndexMatrix(NotNaNIndices), 'off');
if (p < 0.05)
    Comparisons = multcompare(stats, 'display', 'off');
    % only plotting * for anything that is significantly different from
    % undirected
    SigIndices = find((Comparisons(:,end) < 0.05) & (Comparisons(:,2) == length(ConditionDistances)));
    for j = 1:length(SigIndices),
        plot(ConditionDistances(Comparisons(SigIndices(j,1))), 1.15*(Temp(4))/1.2, 'k*', 'MarkerSize', 10, 'LineWidth', 2);
    end
end

% Now calculate correlation to the directed data alone and for the birds that have
% all distances
% Take only rows that have at least 3 not-nan values in the first 3 columns

[Nanrows, Nancols] = find(isnan(PercentChange_DirectedBouts_CV(:,1:3)));
RowsToUse = setdiff(1:1:size(PercentChange_DirectedBouts_CV,1), unique(Nanrows));

DataMatrix = DeltaChange_AllBouts_CV(RowsToUse,1:3);
IndexMatrix = ConditionIndexMatrix(RowsToUse,1:3);

[AllBird_DirectedPearsonR, AllBird_DirectedPearsonP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:));
[AllBird_DirectedSpearmanR, AllBirds_DirectedSpearmanP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:), 'type', 'Spearman', 'rows', 'complete');
text(30, Temp(4)/1.2, ['r = ', num2str(AllBird_DirectedPearsonR), '; p = ', num2str(AllBird_DirectedPearsonP)], 'FontSize', 12, 'FontWeight', 'bold');


% Now to do the same thing for all undirected only bouts

% Normalise as % of undirected and also as delta change relative to
% undirected
for i = 1:size(UnDirectedBouts_CV,1),
    PercentChange_UnDirectedBouts_CV(i,:) = ((UnDirectedBouts_CV(i,:) - UnDirectedBouts_CV(i,end)) * 100)/UnDirectedBouts_CV(i,end);
    DeltaChange_UnDirectedBouts_CV(i,:) = (UnDirectedBouts_CV(i,:) - UnDirectedBouts_CV(i,end));
end

% First for percent change
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [680 289 750 675]);
for i = 1:size(UnDirectedBouts_CV,2),
    AllBout_Bar(i) = bar(ConditionDistances(i), nanmean(PercentChange_UnDirectedBouts_CV(:,i)));
    set(AllBout_Bar(i), 'FaceColor', 'none', 'EdgeColor', 'k', 'BarWidth', 10, 'LineWidth', 2);
    hold on;
    errorbar(ConditionDistances(i), nanmean(PercentChange_UnDirectedBouts_CV(:,i)), nanstd(PercentChange_UnDirectedBouts_CV(:,i))/sqrt(length(find(~isnan(PercentChange_UnDirectedBouts_CV(:,i))))), 'ks-', 'MarkerSize', 10, 'LineWidth', 2);
end
plot(ConditionDistances, nanmean(PercentChange_UnDirectedBouts_CV), 'k', 'LineWidth', 2);

for i = 1:size(UnDirectedBouts_CV,1),
    PercentChangeLine(i) = plot(ConditionDistances, PercentChange_UnDirectedBouts_CV(i,:), [Symbols(mod(i,length(Symbols))+1), Colours(ceil(i/length(Symbols))), '-'], 'LineWidth', 0.5, 'MarkerSize', 5);
end
legend(PercentChangeLine, BirdNames);

clear AxisLabelString;
for j = 1:length(ConditionDistances),
    if (isempty(strfind(IndividualBirds(MaxNumConditionsIndex).Conditions{j}, 'UN')))
        AxisLabelString{j} = num2str(ConditionDistances(j));
    else
        AxisLabelString{j} = 'UN';
    end
end

set(gca, 'XTick', ConditionDistances, 'XTickLabel', AxisLabelString);
xlabel('Distance of female cage (cm)');
ylabel('% change in CV of FF relative to Undirected song');
title(['UnDirected bouts: (n=', num2str(length(BirdNames)), ') birds']);
axis tight;
Temp = axis;
Temp = [(ConditionDistances(1) - 12) (ConditionDistances(end) + 12) 1.02*Temp(3) 1.2*Temp(4)];
axis(Temp);

ConditionIndexMatrix = repmat(1:1:size(PercentChange_UnDirectedBouts_CV,2), size(PercentChange_UnDirectedBouts_CV, 1), 1);
% Omit nans for doing kruskalwallis
%DataMatrix = PercentChange_UnDirectedBouts_CV(RowsToUse,:);
IndexMatrix = ConditionIndexMatrix;
NotNaNIndices = find(~isnan(DataMatrix));
[p, tabl, stats] = kruskalwallis(DataMatrix(NotNaNIndices), IndexMatrix(NotNaNIndices), 'off');
if (p < 0.05)
    Comparisons = multcompare(stats, 'display', 'off');
    % only plotting * for anything that is significantly different from
    % undirected
    SigIndices = find((Comparisons(:,end) < 0.05) & (Comparisons(:,2) == length(ConditionDistances)));
    for j = 1:length(SigIndices),
        plot(ConditionDistances(Comparisons(SigIndices(j,1))), 1.15*(Temp(4))/1.2, 'k*', 'MarkerSize', 10, 'LineWidth', 2);
    end
end

% Now calculate correlation to the directed data alone and for the birds that have
% all distances
% Take only rows that have at least 3 not-nan values in the last 3 female
% presentation columns

[Nanrows, Nancols] = find(isnan(PercentChange_UnDirectedBouts_CV(:,(length(ConditionDistances) - 3):(length(ConditionDistances) - 1))));
RowsToUse = setdiff(1:1:size(PercentChange_UnDirectedBouts_CV,1), unique(Nanrows));

DataMatrix = DeltaChange_AllBouts_CV(RowsToUse,(length(ConditionDistances) - 3):(length(ConditionDistances) - 1));
IndexMatrix = ConditionIndexMatrix(RowsToUse,(length(ConditionDistances) - 3):(length(ConditionDistances) - 1));

[AllBird_DirectedPearsonR, AllBird_DirectedPearsonP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:));
[AllBird_DirectedSpearmanR, AllBirds_DirectedSpearmanP] = corr(ConditionDistances(IndexMatrix(:))', DataMatrix(:), 'type', 'Spearman', 'rows', 'complete');
text(30, Temp(4)/1.2, ['r = ', num2str(AllBird_DirectedPearsonR), '; p = ', num2str(AllBird_DirectedPearsonP)], 'FontSize', 12, 'FontWeight', 'bold');
disp('Finished plotting');