function [] = Harini_PlotFFCV(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% == PCA analysis for all bouts ===========================================
% THe plan is to do a PCA analysis on various features of different bouts.
% Then I can plot directed and undirected bouts on this at different
% distances and then look at the mean position and distances between bouts
% as a measure of how different they are.
% Features that I plan to include are:
% 1) Bout length
% 2) # of INs
% 3) # of syllables
% 4) Time of first motif syllable
% 5) # of motifs
% 6) # of syllables
% 7) Mean FF
% 8) std FF
% 9) Mean motif duration
% 10) std motif duration

% =========================================================================

Colours = 'rgbcmk';
Symbols = '+o<sd^>*pxvh';
Conditions = [{'L0'} {'L1'} {'L2'} {'L3'} {'L4'} {'UN'}];
MinTrialNo = 5;

% PCAFeatures = [{'NumINs'} {'NumSylls'} {'FirstMotifSyllTime'} {'NumMotifs'} {'BoutLength'} {'AverageMotifDur'} {'STDMotifDur'} {'AverageSyllDur'} {'STDSyllDur'} {'AverageSyllEntropy'} {'STDSyllEntropy'} {'AverageSyllMeanFreq'} {'STDSyllMeanFreq'} {'NumTotalSylls'} {'SyllRepetitionRate'}];
% PCAFeatures = [{'FirstMotifSyllTime_GapsLessThan500ms'} {'BoutLength'} {'AverageMotifDur'} {'STDMotifDur'} {'AverageSyllDur'} {'STDSyllDur'} {'AverageSyllEntropy'} {'STDSyllEntropy'} {'AverageSyllMeanFreq'} {'STDSyllMeanFreq'} {'NumTotalSylls'} {'SyllRepetitionRate'}];
PCAFeatures = [{'NumINs'} {'NumMotifs'} {'FirstMotifDur'} {'AverageMotifDur'} {'STDMotifDur'} {'AverageSyllDur'} {'STDSyllDur'} {'AverageSyllEntropy'} {'STDSyllEntropy'} {'AverageSyllMeanFreq'} {'STDSyllMeanFreq'} {'NumTotalSylls'} {'SyllRepetitionRate'}];

for i = 1:length(IndividualBirds),
    for j = 1:length(PCAFeatures),
        ColumnIndex = find(strcmp(PCAFeatures{j}, IndividualBirds(i).BoutStatisticsColumnNames));
        PCAFeatVals{i}(:,j) = IndividualBirds(i).BoutStatistics(:,ColumnIndex);
    end
    
    % Find nan rows
    [NanRows, NanCols] = find(isnan(PCAFeatVals{i}));
    
    % Make all values nans in nan rows
    PCAFeatVals{i}(unique(NanRows),:) = NaN;
    
    % Convert this into a z-score
    PCAFeatVals{i} = (PCAFeatVals{i} - repmat(nanmean(PCAFeatVals{i}), size(PCAFeatVals{i},1),1))./repmat(nanstd(PCAFeatVals{i}), size(PCAFeatVals{i},1),1);
    
    % Now do pca
    [PCA_Coeff, PCA_Score, PCA_Latent, PCA_TSquared, PCA_Explained] = pca(PCAFeatVals{i});

    % First plot explained variance
    figure;
    plot(cumsum(PCA_Explained), 'ko-');
    xlabel('Principal component #');
    ylabel('Cumulative variance explained (%)');
    
    % First plot the co-efficients
    figure;
    hold on;
    for j = 1:6,
        plot(PCA_Coeff(:,j), [Colours(j), 'o-']);
    end
    set(gca, 'XTick', 1:1:size(PCAFeatVals{i}, 2), 'XTickLabel', PCAFeatures, 'XTickLabelRotation', 45)
   
    
    % Next plot the points for different bouts based on the condition
    figure;
    hold on;
    for j = min(IndividualBirds(i).BoutStatistics(:,end)):1:max(IndividualBirds(i).BoutStatistics(:,end)),
        Songs = find((IndividualBirds(i).BoutStatistics(:,end)) == j);
        plot(PCA_Score(Songs, 1), PCA_Score(Songs, 2), [Colours(j), 'o'], 'MarkerFaceColor', Colours(j));
    end
    
    % Next plot the mean and 1 * std dev as confidence ellipses for the data
    figure;
    hold on;
    for j = min(IndividualBirds(i).BoutStatistics(:,end)):1:max(IndividualBirds(i).BoutStatistics(:,end)),
        Songs = find(((IndividualBirds(i).BoutStatistics(:,end)) == j) & (~isnan(PCA_Score(:,1))));
        PlotConfidenceEllipse(PCA_Score(Songs,1:2), Colours(j), 1)
        plot(nanmean(PCA_Score(Songs, 1)), nanmean(PCA_Score(Songs, 2)), [Colours(j), '+'], 'MarkerSize', 10, 'LineWidth', 2);
    end
end

disp('Finished PCA analysis on Bout features');