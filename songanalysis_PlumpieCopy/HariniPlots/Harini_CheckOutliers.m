function [] = Harini_CheckOutliers(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

Fid = fopen('OutliersLog.txt', 'w');
for i = 1:length(IndividualBirds),
    fprintf(Fid, '%s\n', BirdNames{i});
    
    Templates = load(IndividualBirds(i).SortedBirdParameters(1).TemplateFile);
    for j = 1:length(Templates.SyllableTemplates),
        SyllablesWithTemplates(j) = Templates.SyllableTemplates{j}{1}.MotifTemplate(1).Label;
    end
    
    for j = 1:length(SyllablesWithTemplates),
        % First find matches
        Matches = find((char(IndividualBirds(i).AllSyllableData(:,1)) == SyllablesWithTemplates(j)));
        
        % Next find outliers based on template match values
        TemplateMatchOutlierThreshold(1) = (prctile(IndividualBirds(i).AllSyllableTemplateMatchValues(Matches), 25) - 2*iqr(IndividualBirds(i).AllSyllableTemplateMatchValues(Matches)));
        TemplateMatchOutlierThreshold(2) = (prctile(IndividualBirds(i).AllSyllableTemplateMatchValues(Matches), 75) + 2*iqr(IndividualBirds(i).AllSyllableTemplateMatchValues(Matches)));
        fprintf(Fid, 'Syll %c\tTemplate match Outlier thresholds = %f and %f\n', SyllablesWithTemplates(j), TemplateMatchOutlierThreshold(1), TemplateMatchOutlierThreshold(2));
        
        Outliers = find((IndividualBirds(i).AllSyllableTemplateMatchValues(Matches) < TemplateMatchOutlierThreshold(1)) | (IndividualBirds(i).AllSyllableTemplateMatchValues(Matches) > TemplateMatchOutlierThreshold(2)));
        Outliers = Matches(Outliers);
        
        TempMatchOutliers = Outliers;
        TempMatchOutlierValues = IndividualBirds(i).AllSyllableTemplateMatchValues(TempMatchOutliers);
        
        FeatureCols = [1 3 4 6 7 8]; % exclude amplitude and amplitude modulation
        
        % Now find outliers based on mahalanobis distance
        [NanRows, NanCols] = find(isnan(IndividualBirds(i).AllSyllableFeatValues(Matches,FeatureCols)));
        NanRows = unique(NanRows);
        NanRows = Matches(NanRows);
        NonNanRows = setdiff(Matches, NanRows);
        if (length(NonNanRows) > 8)
            Distances = pdist2(IndividualBirds(i).AllSyllableFeatValues(NonNanRows,FeatureCols), mean(IndividualBirds(i).AllSyllableFeatValues(NonNanRows,FeatureCols)), 'mahalanobis', cov(IndividualBirds(i).AllSyllableFeatValues(NonNanRows,FeatureCols)));
        
            DistanceOutlierThreshold(1) = (prctile(Distances, 25) - 2*iqr(Distances));
            DistanceOutlierThreshold(2) = (prctile(Distances, 75) + 2*iqr(Distances));
        
            fprintf(Fid, 'Syll %c\tMahalanobis Distance Outlier thresholds = %f and %f\n', SyllablesWithTemplates(j), DistanceOutlierThreshold(1), DistanceOutlierThreshold(2));
        
            Outliers = find((Distances < DistanceOutlierThreshold(1)) | (Distances > DistanceOutlierThreshold(2)));
            DistanceOutlierIndices = NonNanRows(Outliers);
            OutlierDistances = Distances(Outliers);
            
            AllOutliers = [TempMatchOutliers(:); DistanceOutlierIndices(:)];
            AllOutlierValues = [TempMatchOutlierValues(:); OutlierDistances(:)];
            AllOutlierIdentity = [ones(size(TempMatchOutliers(:))); ones(size(DistanceOutlierIndices(:)))*2];
            OutlierIdentityString = {'Template match' 'Mahalanobis distance'};
            
            [UniqueOutliers, UniqueOutlierIndices, IC] = unique(AllOutliers);
            fprintf(Fid, '>>> # of outliers = %i/%i; %f%%\n', length(UniqueOutlierIndices), length(Matches), 100*length(UniqueOutlierIndices)/length(Matches));
            for k = 1:length(UniqueOutlierIndices),
                fprintf(Fid, 'Outlier #%i: File Name. %s: Syll onset time %fms; Outlier type %s; Match value = %f\n', k, IndividualBirds(i).SongFileNames{IndividualBirds(i).AllSyllableData(AllOutliers(UniqueOutlierIndices(k)),2)}, IndividualBirds(i).AllSyllableData(AllOutliers(UniqueOutlierIndices(k)),4), OutlierIdentityString{AllOutlierIdentity(UniqueOutlierIndices(k))}, AllOutlierValues(UniqueOutlierIndices(k)));
            end
        else
            fprintf(Fid, '>>> # of outliers = %i/%i; %f%%\n', length(Outliers), length(Matches), 100*length(Outliers)/length(Matches));
            for k = 1:length(Outliers),
                fprintf(Fid, 'Outlier #%i: File Name: %s: Syll onset time %fms; Template Match; Match value = %f\n', k, IndividualBirds(i).SongFileNames{IndividualBirds(i).AllSyllableData(Outliers(k),2)}, IndividualBirds(i).AllSyllableData(Outliers(k),4), Distances(Outliers(k))); 
            end
        end
    end
    fprintf(Fid, '\n');
end

fclose(Fid);

disp('Finished checking outliers');
