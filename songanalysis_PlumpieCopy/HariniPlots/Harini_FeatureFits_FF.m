function [L0Dir_Undir_PValue] = Harini_FeatureFits_FF(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, Distances, BirdNames)

FinalFigureDir = '/home/raghav/StudentRelated/Harini/Manuscript/results';
DistanceCode = {'L0' 'L1' 'L2' 'L3' 'L4' 'UN'};

% We want to do the analysis per bird and for the group as a whole. The
% things that we will analyze are (1) # of INs at the beginning based on
% the 500ms criterion, (2) # of complete motifs in the bout, (3) all and
% first motif duration, (4), CV of FF and (5) amplitude of individual
% motif syllables only for HF mic birds

MinTrialNo = 7;
% Now to do the log amplitude calculated similar to Kao et al.
for i = 1:length(DirBout_Stats)
    disp(['Bird #', num2str(i), ': ', BirdNames{i}]);
    % GroupData has 3 columns, the first one has the data corresponding to
    % the variable that that I'm looking at like TotalINNumber, etc. 2nd
    % column has the index of the distance at which the recordings were
    % done (0 for L0, 1 for L1 ... 6 for UN) and 3rd column has 1 for Dir,
    % 2 for Dir/Undir - DUN - and 3 for Undir
    if (~isfield(UndirBout_Stats{i}{end}, 'FF_SyllLabel'))
        continue;
    end
    figure(i);
    FigPanel(i).p = panel();
    FigPanel(i).p.pack(1, length(UndirBout_Stats{i}{end}.FF_SyllLabel));
    
    for SyllIndex = 1:length(UndirBout_Stats{i}{end}.FF_SyllLabel),
                
        FFGroupData{i}{SyllIndex} = [];
        SyllDurationGroupData{i}{SyllIndex} = [];
        for j = 1:length(DirBout_Stats{i}),
            if (~isempty(DirBout_Stats{i}{j}))
                FFGroupData{i}{SyllIndex} = [FFGroupData{i}{SyllIndex}; [DirBout_Stats{i}{j}.FF{SyllIndex}(:) ones(size(DirBout_Stats{i}{j}.FF{SyllIndex}(:)))*j ones(size(DirBout_Stats{i}{j}.FF{SyllIndex}(:)))*1]];
                SyllDurationGroupData{i}{SyllIndex} = [SyllDurationGroupData{i}{SyllIndex}; [DirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:) ones(size(DirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:)))*j ones(size(DirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:)))*1]];
            end
        end
    
        for j = 1:length(DirUndirBout_Stats{i}),
            if (~isempty(DirUndirBout_Stats{i}{j}))
                FFGroupData{i}{SyllIndex} = [FFGroupData{i}{SyllIndex}; [DirUndirBout_Stats{i}{j}.FF{SyllIndex}(:) ones(size(DirUndirBout_Stats{i}{j}.FF{SyllIndex}(:)))*j ones(size(DirUndirBout_Stats{i}{j}.FF{SyllIndex}(:)))*2]];
                SyllDurationGroupData{i}{SyllIndex} = [SyllDurationGroupData{i}{SyllIndex}; [DirUndirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:) ones(size(DirUndirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:)))*j ones(size(DirUndirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:)))*2]];
            end
        end
        
        for j = 1:length(UndirBout_Stats{i}),
            if (~isempty(UndirBout_Stats{i}{j}))
                FFGroupData{i}{SyllIndex} = [FFGroupData{i}{SyllIndex}; [UndirBout_Stats{i}{j}.FF{SyllIndex}(:) ones(size(UndirBout_Stats{i}{j}.FF{SyllIndex}(:)))*j ones(size(UndirBout_Stats{i}{j}.FF{SyllIndex}(:)))*3]];
                SyllDurationGroupData{i}{SyllIndex} = [SyllDurationGroupData{i}{SyllIndex}; [UndirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:) ones(size(UndirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:)))*j ones(size(UndirBout_Stats{i}{j}.SyllDuration{SyllIndex}(:)))*3]];
            end
        end

        % Remove NaN values
        NaNIndices = find(isnan(FFGroupData{i}{SyllIndex}(:,1)));
        if (~isempty(NaNIndices))
            FFGroupData{i}{SyllIndex}(NaNIndices,:) = [];
            SyllDurationGroupData{i}{SyllIndex}(NaNIndices,:) = [];
        end
        
        % Remove outliers based on the criterion that data is either 3*IQR 
        % above 75th percentile or below 25th percentile for syllable
        % duration and then for amplitude
        OutlierThreshold = [(prctile(SyllDurationGroupData{i}{SyllIndex}(:,1), 25) - 3*iqr(SyllDurationGroupData{i}{SyllIndex}(:,1))) (prctile(SyllDurationGroupData{i}{SyllIndex}(:,1), 75) + 3*iqr(SyllDurationGroupData{i}{SyllIndex}(:,1)))];
        OutlierIndices = find((SyllDurationGroupData{i}{SyllIndex}(:,1) < OutlierThreshold(1)) | (SyllDurationGroupData{i}{SyllIndex}(:,1) > OutlierThreshold(2)));
        if (~isempty(OutlierIndices))
            SyllDurationGroupData{i}{SyllIndex}(OutlierIndices, :) = [];
            FFGroupData{i}{SyllIndex}(OutlierIndices, :) = [];
            disp(['Removed ', num2str(length(OutlierIndices)), ' from syllable ', UndirBout_Stats{i}{end}.FF_SyllLabel(SyllIndex), ' based on syll duration:', BirdNames{i}]);
        end
    
        % Remove outliers based on the criterion that data is either 3*IQR 
        % above 75th percentile or below 25th percentile for syllable
        % duration and then for amplitude
        OutlierThreshold = [(prctile(FFGroupData{i}{SyllIndex}(:,1), 25) - 3*iqr(FFGroupData{i}{SyllIndex}(:,1))) (prctile(FFGroupData{i}{SyllIndex}(:,1), 75) + 3*iqr(FFGroupData{i}{SyllIndex}(:,1)))];
        OutlierIndices = find((FFGroupData{i}{SyllIndex}(:,1) < OutlierThreshold(1)) | (FFGroupData{i}{SyllIndex}(:,1) > OutlierThreshold(2)));
        if (~isempty(OutlierIndices))
            SyllDurationGroupData{i}{SyllIndex}(OutlierIndices, :) = [];
            FFGroupData{i}{SyllIndex}(OutlierIndices, :) = [];
            disp(['Removed ', num2str(length(OutlierIndices)), ' from syllable ', UndirBout_Stats{i}{end}.FF_SyllLabel(SyllIndex), ' based on syll amplitude:', BirdNames{i}]);
        end
        
        % Now find and remove any distance that doesn't have 3 bouts of 
        % directed song
        for j = 1:max(FFGroupData{i}{SyllIndex}(:,end-1)),
            DirIndicesForDistance = find((FFGroupData{i}{SyllIndex}(:,end-1) == j) & (FFGroupData{i}{SyllIndex}(:,end) == 1));
            if ((length(DirIndicesForDistance) > 0) & (length(DirIndicesForDistance) < MinTrialNo))
                FFGroupData{i}{SyllIndex}(DirIndicesForDistance,:) = [];
                disp(['Removed distance ', DistanceCode{j}, ' for ', BirdNames{i}, ' as there were only ', num2str(length(DirIndicesForDistance)), ' syllables for syllable ', UndirBout_Stats{i}{end}.FF_SyllLabel(SyllIndex), ' - log ampitude']);
            end
        end
    
        DirIndices = find(FFGroupData{i}{SyllIndex}(:,end) == 1);
        UnDirIndices = find((FFGroupData{i}{SyllIndex}(:,end-1) == 6) & (FFGroupData{i}{SyllIndex}(:,end) == 3));
        % First plot the group data
        FigPanel(i).p(1,SyllIndex).select();
        hold on;

        % Now plot the mean, SEM and the best fit line with robust regression
        CV = [];
        for Location = 1:max(FFGroupData{i}{SyllIndex}(DirIndices,end-1)),
            LocIndices = intersect(DirIndices, find(FFGroupData{i}{SyllIndex}(:,end-1) == Location));
            if (~isempty(LocIndices))
                plot(Distances(Location), std(FFGroupData{i}{SyllIndex}(LocIndices, 1))/mean(FFGroupData{i}{SyllIndex}(LocIndices, 1)), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
                CV(end+1,:) = [Distances(Location), std(FFGroupData{i}{SyllIndex}(LocIndices, 1))/mean(FFGroupData{i}{SyllIndex}(LocIndices, 1))];
            end
        end
        plot(Distances(Location)+30, std(FFGroupData{i}{SyllIndex}(UnDirIndices, 1))/mean(FFGroupData{i}{SyllIndex}(UnDirIndices, 1)), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
        
        % Now best fit line
        RobustCoeffs = robustfit(CV(:,1)', CV(:,2)');

        [FixedX_BootStrapCoeffs, Random_BootStrapCoeffs, RandomProb, Rsq] = BootstrapTestRegressionSignificance(RobustCoeffs, CV(:,1), CV(:,2), 'linear');

        if (RobustCoeffs(2) < 0)
            if ((FixedX_BootStrapCoeffs(1,2) < 0) && (FixedX_BootStrapCoeffs(3,2) < 0))
                PlotLineColor = 'b';
            else
                PlotLineColor = 'k';
            end
        else
            if (RobustCoeffs(2) > 0)
                if ((FixedX_BootStrapCoeffs(1,2) > 0) && (FixedX_BootStrapCoeffs(3,2) > 0))
                    PlotLineColor = 'r';
                else
                    PlotLineColor = 'k';
                end
            else
                PlotLineColor = 'k';
            end
        end
        plot(unique(Distances(FFGroupData{i}{SyllIndex}(DirIndices,end-1))), polyval(flipud(RobustCoeffs), CV(:,1)), PlotLineColor, 'LineWidth', 1.5);
        patch([CV(:,1) fliplr(CV(:,1))], [polyval(fliplr(FixedX_BootStrapCoeffs(1,:)), CV(:,1)) fliplr(polyval(fliplr(FixedX_BootStrapCoeffs(3,:)), CV(:,1)))], PlotLineColor, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        
        axis tight;
        Temp = axis;
        Temp = [-5 Temp(2)+5 0.99*Temp(3) 1.01*Temp(4)];
        text(Temp(2)-23, Temp(3), '//', 'FontSize', 16);
        plot(Temp(1:2), ones(1,2)*RobustCoeffs(1), 'k--');
        axis(Temp);
        ylabel('CV of FF')
        xlabel('Distance (cm)');

    end
end


disp('Done with stats and plots');
