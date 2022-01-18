function [L0Dir_Undir_PValue] = Harini_AllFeatureStats_Plots(DirBout_Stats, DirUndirBout_Stats, UndirBout_Stats, StructString, MinTrialNo, Distances, LabelString, BirdNames)

% Now to do the stats
% Now, we want to normalize all of the features across all the data for
% each bird. Then plot the normalized changes for directed songs for all
% features for a given bird in the same plot.
% Normalization is just a z-score across all data for that bird and that
% feature

FeaturesToConsider = {'TotalINNumber_500ms' 'CompleteMotifNumber' 'FirstMotifDuration'};

Colours = 'rbkcmg';

for i = 1:length(DirBout_Stats)
    disp(['Bird #', num2str(i), ': ', BirdNames{i}]);
    GroupData{i} = [];
    % GroupData has 3 columns, the first one has the data corresponding to
    % the variable that that I'm looking at like TotalINNumber, etc. 2nd
    % column has the index of the distance at which the recordings were
    % done (0 for L0, 1 for L1 ... 6 for UN) and 3rd column has 1 for Dir,
    % 2 for Dir/Undir - DUN - and 3 for Undir
    
    for j = 1:length(DirBout_Stats{i}),
        if (~isempty(DirBout_Stats{i}{j}))
            TempGroupData = [];
            for k = 1:length(FeaturesToConsider),
                if (k == length(FeaturesToConsider))
                    TempGroupData = [TempGroupData [eval(['DirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)']) ones(size(eval(['DirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])))*j ones(size(eval(['DirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])))*1]];
                else
                    TempGroupData = [TempGroupData eval(['DirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])];
                end
            end
            GroupData{i} = [GroupData{i}; TempGroupData];
        end
    end
    
    for j = 1:length(DirUndirBout_Stats{i}),
        if (~isempty(DirUndirBout_Stats{i}{j}))
            TempGroupData = [];
            for k = 1:length(FeaturesToConsider),
                if (k == length(FeaturesToConsider))
                    TempGroupData = [TempGroupData [eval(['DirUndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)']) ones(size(eval(['DirUndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])))*j ones(size(eval(['DirUndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])))*2]];
                else
                    TempGroupData = [TempGroupData eval(['DirUndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])];
                end
            end
            GroupData{i} = [GroupData{i}; TempGroupData];
        end
    end
    
    for j = 1:length(UndirBout_Stats{i}),
        if (~isempty(UndirBout_Stats{i}{j}))
            TempGroupData = [];
            for k = 1:length(FeaturesToConsider),
                if (k == length(FeaturesToConsider))
                    TempGroupData = [TempGroupData [eval(['UndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)']) ones(size(eval(['UndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])))*j ones(size(eval(['UndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])))*3]];
                else
                    TempGroupData = [TempGroupData eval(['UndirBout_Stats{i}{j}.', FeaturesToConsider{k}, '(:)'])];
                end
            end
            GroupData{i} = [GroupData{i}; TempGroupData];
        end
    end

    % Remove NaN values
    NaNIndices = find(isnan(GroupData{i}(:,1)));
    if (~isempty(NaNIndices))
        GroupData{i}(NaNIndices,:) = [];
    end
    
    % GroupData{i}(:,1:length(FeaturesToConsider)-1) = log10(GroupData{i}(:,1:length(FeaturesToConsider)-1));
    % GroupData{i}(:,1:length(FeaturesToConsider)) = zscore(GroupData{i}(:,1:length(FeaturesToConsider)));
    [PCACoeff{i}, PCAScore{i}, PCALatent{i}] = pca((GroupData{i}(:,1:length(FeaturesToConsider))));
    
    % Now plot the distributions and means for all birds across all
    % distances
end

% Now put all birds together and calculateo pca on all birds so that I can
% keep the same coeffs for all birds
AllBirdData = [];
for i = 1:length(GroupData),
    AllBirdData = [AllBirdData; GroupData{i}];
end
[AllBirdPCACoeff, AllBirdPCAScore, AllBirdPCALatent] = pca((AllBirdData(:,1:length(FeaturesToConsider))));


% Now to plot the mean and std across distances for D, D/UN and UN for all
% birds together as percent change relative to undirected.
% This is done only for birds where the means are already different between
% L0 directed and undirected as measured by the earlier p-value

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
    L0DirSongIndices = find((GroupData{i}(:,length(FeaturesToConsider)+2) == 1) & (GroupData{i}(:,length(FeaturesToConsider)+1) == 1));
    UnDirUnDirSongIndices = find((GroupData{i}(:,length(FeaturesToConsider)+2) == 3) & (GroupData{i}(:,length(FeaturesToConsider)+1) == 6));
    
    % GroupData{i}(:,1:length(FeaturesToConsider)) = zscore(GroupData{i}(:,1:length(FeaturesToConsider)));
    
    % Recalculate PCA scores based on common axis across all birds
    PCAScore{i} = GroupData{i}(:,1:length(FeaturesToConsider)) * AllBirdPCACoeff;
    
    % First for directed songs
    for j = 1:6,
        DirSongIndices = find((GroupData{i}(:,length(FeaturesToConsider)+2) == 1) & (GroupData{i}(:,length(FeaturesToConsider)+1) == j));
        if (length(DirSongIndices) >= MinTrialNo)
            for k = 1:length(FeaturesToConsider),
                DistanceDirMeans{k}(Index, j) = mean(GroupData{i}(DirSongIndices,k));
                DistanceDirSEMs{k}(Index, j) = std(GroupData{i}(DirSongIndices,k))/sqrt(length(DirSongIndices));
                DistanceDirCIs{k}{Index}(j,:) = bootci(10000,@mean, GroupData{i}(DirSongIndices,k));
            end
            if (length(L0DirSongIndices) >= MinTrialNo)
                %MahalDistances_L0Dir{Index}{j} = pdist2(GroupData{i}(DirSongIndices, 1:length(FeaturesToConsider)), mean(GroupData{i}(L0DirSongIndices, 1:length(FeaturesToConsider))), 'mahal', cov(GroupData{i}(L0DirSongIndices, 1:length(FeaturesToConsider))));
                EucDistances_L0Dir{Index}{j} = pdist2(PCAScore{i}(DirSongIndices, 1:length(FeaturesToConsider)), mean(PCAScore{i}(L0DirSongIndices, 1:length(FeaturesToConsider))));
                
                %MeanMahalDistances_L0Dir(Index,j) = mean(MahalDistances_L0Dir{Index}{j});
                MeanEucDistances_L0Dir(Index,j) = mean(EucDistances_L0Dir{Index}{j});
                
                %CI = bootci(10000, @mean, MahalDistances_L0Dir{Index}{j});
                %MahalLowerCI_L0Dir(Index,j) = CI(1) - mean(MahalDistances_L0Dir{Index}{j});
                %MahalUpperCI_L0Dir(Index,j) = CI(2) - mean(MahalDistances_L0Dir{Index}{j});
                
                CI = bootci(10000, @mean, EucDistances_L0Dir{Index}{j});
                EucLowerCI_L0Dir(Index,j) = CI(1) - mean(EucDistances_L0Dir{Index}{j});
                EucUpperCI_L0Dir(Index,j) = CI(2) - mean(EucDistances_L0Dir{Index}{j});
            else
%                 MahalDistances_L0Dir{Index}{j} = NaN;
%                 MeanMahalDistances_L0Dir(Index,j) = NaN;
%                 MahalLowerCI_L0Dir(Index,j) = NaN;
%                 MahalUpperCI_L0Dir(Index,j) = NaN;

                EucDistances_L0Dir{Index}{j} = NaN;
                MeanEucDistances_L0Dir(Index,j) = NaN;
                EucLowerCI_L0Dir(Index,j) = NaN;
                EucUpperCI_L0Dir(Index,j) = NaN;
            end
            
            if (length(UnDirUnDirSongIndices) >= MinTrialNo)
                %MahalDistances_UnDirUnDir{Index}{j} = pdist2(GroupData{i}(DirSongIndices, 1:length(FeaturesToConsider)), mean(GroupData{i}(UnDirUnDirSongIndices, 1:length(FeaturesToConsider))), 'mahal', cov(GroupData{i}(UnDirUnDirSongIndices, 1:length(FeaturesToConsider))));
                EucDistances_UnDirUnDir{Index}{j} = pdist2(PCAScore{i}(DirSongIndices, 1:length(FeaturesToConsider)), mean(PCAScore{i}(UnDirUnDirSongIndices, 1:length(FeaturesToConsider))));
                
                %MeanMahalDistances_UnDirUnDir(Index,j) = mean(MahalDistances_UnDirUnDir{Index}{j});
                MeanEucDistances_UnDirUnDir(Index,j) = mean(EucDistances_UnDirUnDir{Index}{j});
                
%                 CI = bootci(10000, @mean, MahalDistances_UnDirUnDir{Index}{j});
%                 MahalLowerCI_UnDirUnDir(Index,j) = CI(1) - mean(MahalDistances_UnDirUnDir{Index}{j});
%                 MahalUpperCI_UnDirUnDir(Index,j) = CI(2) - mean(MahalDistances_UnDirUnDir{Index}{j});
                
                CI = bootci(10000, @mean, EucDistances_UnDirUnDir{Index}{j});
                EucLowerCI_UnDirUnDir(Index,j) = CI(1) - mean(EucDistances_UnDirUnDir{Index}{j});
                EucUpperCI_UnDirUnDir(Index,j) = CI(2) - mean(EucDistances_UnDirUnDir{Index}{j});
            else
%                 MahalDistances_UnDirUnDir{Index}{j} = NaN;
%                 MeanMahalDistances_UnDirUnDir(Index,j) = NaN;
%                 MahalLowerCI_UnDirUnDir(Index,j) = NaN;
%                 MahalUpperCI_UnDirUnDir(Index,j) = NaN;    
                
                EucDistances_UnDirUnDir{Index}{j} = NaN;
                MeanEucDistances_UnDirUnDir(Index,j) = NaN;
                EucLowerCI_UnDirUnDir(Index,j) = NaN;
                EucUpperCI_UnDirUnDir(Index,j) = NaN;    
            end
        else
            for k = 1:length(FeaturesToConsider),
                DistanceDirMeans{k}(Index, j) = NaN;
                DistanceDirSEMs{k}(Index,j) = NaN;
                DistanceDirCIs{k}{Index}(j,:) = ones(1,2)*NaN;
            end
%             MahalDistances_L0Dir{Index}{j} = NaN;
%             MeanMahalDistances_L0Dir(Index,j) = NaN;
%             MahalLowerCI_L0Dir(Index,j) = NaN;
%             MahalUpperCI_L0Dir(Index,j) = NaN;
            
%             MahalDistances_UnDirUnDir{Index}{j} = NaN;
%             MeanMahalDistances_UnDirUnDir(Index,j) = NaN;
%             MahalLowerCI_UnDirUnDir(Index,j) = NaN;
%             MahalUpperCI_UnDirUnDir(Index,j) = NaN;
            
            EucDistances_L0Dir{Index}{j} = NaN;
            MeanEucDistances_L0Dir(Index,j) = NaN;
            EucLowerCI_L0Dir(Index,j) = NaN;
            EucUpperCI_L0Dir(Index,j) = NaN;
            
            EucDistances_UnDirUnDir{Index}{j} = NaN;
            MeanEucDistances_UnDirUnDir(Index,j) = NaN;
            EucLowerCI_UnDirUnDir(Index,j) = NaN;
            EucUpperCI_UnDirUnDir(Index,j) = NaN;
        end
        
        DirUnDirSongIndices = find((GroupData{i}(:,length(FeaturesToConsider)+2) == 2) & (GroupData{i}(:,length(FeaturesToConsider)+1) == j));
        if (length(DirUnDirSongIndices) >= MinTrialNo)
            for k = 1:length(FeaturesToConsider),
                DistanceDirUnDirMeans{k}(Index, j) = mean(GroupData{i}(DirUnDirSongIndices,k));
                DistanceDirUnDirSEMs{k}(Index, j) = std(GroupData{i}(DirUnDirSongIndices,k))/sqrt(length(DirUnDirSongIndices));
                DistanceDirUnDirCIs{k}{Index}(j,:) = bootci(10000, @mean, GroupData{i}(DirUnDirSongIndices,k));
            end
        else
            for k = 1:length(FeaturesToConsider),
                DistanceDirUnDirMeans{k}(Index, j) = NaN;
                DistanceDirUnDirSEMs{k}(Index, j) = NaN;
                DistanceDirUnDirCIs{k}{Index}(j,:) = ones(1,2)*NaN;
            end
        end
        
        UnDirSongIndices = find((GroupData{i}(:,length(FeaturesToConsider)+2) == 3) & (GroupData{i}(:,length(FeaturesToConsider)+1) == j));
        if (length(UnDirSongIndices) >= MinTrialNo)
            for k = 1:length(FeaturesToConsider),
                DistanceUnDirMeans{k}(Index, j) = mean(GroupData{i}(UnDirSongIndices,k));
                DistanceUnDirSEMs{k}(Index, j) = std(GroupData{i}(UnDirSongIndices,k))/sqrt(length(UnDirSongIndices));
                DistanceUnDirCIs{k}{Index}(j,:) = bootci(10000, @mean, GroupData{i}(UnDirSongIndices,k));
            end
        else
            for k = 1:length(FeaturesToConsider),
                DistanceUnDirMeans{k}(Index, j) = NaN;
                DistanceUnDirSEMs{k}(Index, j) = NaN;
                DistanceUnDirCIs{k}{Index}(j,:) = ones(1,2)*NaN;
            end
        end
    end
    
    figure(DirUnDirComparisonFig);
    subplot(3, ceil(length(DirBout_Stats)/3),Index);
    hold on;
    for k = 1:length(FeaturesToConsider),
        plot(Distances, DistanceDirMeans{k}(Index,:), [Colours(k), 's-'], 'LineWidth', 1);
        plot(repmat(Distances(:), 1, 2)', DistanceDirCIs{k}{Index}', Colours(k))
      
        % plot(Distances, DistanceUnDirMeans{k}(Index,:), 'bs-', 'LineWidth', 1);
        %plot(repmat(Distances(:), 1, 2)', DistanceUnDirCIs{Index}', 'b');
        plot(repmat(Distances(:),1,2), repmat(DistanceUnDirCIs{k}{Index}(end,:), length(Distances), 1), [Colours(k), '--']);
    end
    
    xlabel('Distance from female (cm)');
    axis tight;
    Temp = axis;
    Temp = [-0.5 250 0.95*Temp(3) 1.05*Temp(4)];
    axis(Temp);
    ylabel(LabelString);
    title(BirdNames{i}, 'Color', 'r');
    % legend(FeaturesToConsider);
end



disp('Done with stats and plots');
