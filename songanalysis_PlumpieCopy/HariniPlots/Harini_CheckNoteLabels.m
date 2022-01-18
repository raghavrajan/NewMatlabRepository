function [] = Harini_CheckNoteLabels(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% ============== Function to check all note labels ========================
% A general problem associated with labelling is we don't have a measure of
% how different individual clusters are from each other. The first thing
% to do for this would be to plot the variance of each cluster measured as
% the square root of the mahalanobis distance of each point to the centroid
% of that cluster and comparing this with inter-cluster distances.

CutOffPercentage = 1; % this represents the cutoff percentage for relative proportion of each syllable type. Only if the syllable type proportion is >= cut-off, will it be considered

for i = 1:length(IndividualBirds),
    disp(['Bird #', num2str(i)]);
    AllSyllableData = IndividualBirds(i).AllSyllableData;
    AllSyllableFeatValues = IndividualBirds(i).AllSyllableFeatValues;
    
    TotalNumSylls = size(AllSyllableData, 1);
    
    % First remove NaN values from syll feat values and display the % of
    % syllables that have been removed
    [NanRows, NanCols] = find(isnan(AllSyllableFeatValues));
    SyllsToRemove = unique(NanRows);
    disp(['Removed ', num2str(length(SyllsToRemove)), ' syllables out of ', num2str(TotalNumSylls), '(', num2str(100 * length(SyllsToRemove) / TotalNumSylls), '%)']);
    
    AllSyllableData(SyllsToRemove,:) = [];
    AllSyllableFeatValues(SyllsToRemove,:) = [];

    UniqueSyllLabels = unique(char(AllSyllableData(:,1)));
    TotalNumSylls = size(AllSyllableData, 1);
    
    disp(['Total # of syllables = ', num2str(length(UniqueSyllLabels))]);
    
    InterClustDistances{i} = ones(length(UniqueSyllLabels))*NaN;
    
    for j = 1:length(UniqueSyllLabels),
        Sylls = find(char(AllSyllableData(:,1)) == UniqueSyllLabels(j));
        NumSylls = length(Sylls);
        if ((100*NumSylls/TotalNumSylls) >= CutOffPercentage)
            Centroid{i}(j,:) = mean(AllSyllableFeatValues(Sylls,1:8));
            Variance(i,j) = sqrt(sum(pdist2(AllSyllableFeatValues(Sylls,1:8), Centroid{i}(j,:), 'mahalanobis')));
            IndividualDistances{i}{j} = pdist2(AllSyllableFeatValues(Sylls,1:8), Centroid{i}(j,:), 'mahalanobis');
            for k = 1:length(UniqueSyllLabels),
                if (k == j)
                    continue;
                else
                    Sylls2 = find(char(AllSyllableData(:,1)) == UniqueSyllLabels(k));
                    NumSylls2 = length(Sylls2);
                    if ((100*NumSylls2/TotalNumSylls) >= CutOffPercentage)
                        Centroid2 = mean(AllSyllableFeatValues(Sylls2,1:8));
                        InterClustDistances{i}(j,k) = mean(pdist2(AllSyllableFeatValues(Sylls,1:8), Centroid2, 'mahalanobis'));
                    end
                end
            end
        end
        disp(['Syll ', UniqueSyllLabels(j), ': proportion of total syllables = ', num2str(NumSylls * 100 / TotalNumSylls), '% (', num2str(NumSylls), ' / ', num2str(TotalNumSylls), ')']);
    end
end
disp('Finished assessing note labels');