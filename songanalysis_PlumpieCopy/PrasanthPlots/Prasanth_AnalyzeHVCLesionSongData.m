function [] = Prasanth_AnalyzeHVCLesionSongData(BirdDetailsTextFile, MinTrialNo)

% First get details from the CSV text file
disp('Getting header data from CSV file ...');
[HeaderLine, BirdDetails] = LSINA_GetDetailsFromCSVFile(BirdDetailsTextFile);

% Now parse all the lines into the appropriate variables based on the
% header line
disp('Getting data from CSV file ...');
[BirdParameters] = LSINA_ParseHeaderBirdData(HeaderLine, BirdDetails);

% Now for each of the birds, load up all the filenames
disp('Loading up filenames ...');
for i = 1:length(BirdParameters),
    fprintf('%d >> ', i);
    [BirdParameters(i).SongFileNames] = LSINA_GetDataFileNames(BirdParameters(i));
end
fprintf('\n');

% Now load up the note files and the length of each file
disp('Loading up note data ...');
for i = 1:length(BirdParameters),
    fprintf('%d >> ', i);
    [BirdParameters(i).NoteInfo, BirdParameters(i).FileLen] = LSINA_LoadNoteFileInfo(BirdParameters(i));
end
fprintf('\n');

% Now the first thing to do would be to check if it is continuous data or
% not. If it is continuous data, then check for consecutive Capital letter
% syllables that would correspond to the same syllable split over two
% conseecutive files. These have to be merged.
% Do this for continuous data and then put together one long list of
% syllables and their corresponding file #s, onsets and offsets.

disp('Putting together list of syllables ...');
for i = 1:length(BirdParameters),
    fprintf('%d >> ', i);
    [BirdParameters(i).SyllableData] = GetSyllableListInfo(BirdParameters(i));
end
fprintf('\n');

% Now calculate SAP features for all the files
disp('Calculating SAP features ...');
for i = 1:length(BirdParameters),
    disp(['Calculating SAP features for ', BirdParameters(i).BirdName, '-', BirdParameters(i).DataLabel, ' ...']);
    if (exist([BirdParameters(i).BirdName, '.', BirdParameters(i).DataLabel, '.SAPFeats.mat'], 'file'))
        load([BirdParameters(i).BirdName, '.', BirdParameters(i).DataLabel, '.SAPFeats.mat']);
        BirdParameters(i).SAPFeatsMatrix = SAPFeatsMatrix;
        BirdParameters(i).SAPFeat_FieldNames = SAPFeat_FieldNames;
    else
        [BirdParameters(i).SAPFeatsMatrix, BirdParameters(i).SAPFeat_FieldNames] = CalcSAPFeats(BirdParameters(i));
        SAPFeatsMatrix = BirdParameters(i).SAPFeatsMatrix;
        SAPFeat_FieldNames = BirdParameters(i).SAPFeat_FieldNames;
        save([BirdParameters(i).BirdName, '.', BirdParameters(i).DataLabel, '.SAPFeats.mat'], 'SAPFeatsMatrix', 'SAPFeat_FieldNames');
    end
end

% Now split up the files into bouts based on inter-bout interval that is
% also specified in the .csv file
disp('Identifying bouts ...');
for i = 1:length(BirdParameters),
    fprintf('%d >> ', i);
    [BirdParameters(i).Bouts] = IdentifyBouts(BirdParameters(i));
end
fprintf('\n');

% Now to plot the duration of syllables as a function of the position in a
% bout 
disp('Plotting syllable duration vs. bout position ...');
for i = 1:length(BirdParameters),
    fprintf('%d >> ', i);
    IndividualBoutSyllDurations{i} = [];
    IndividualBoutSyllAmplitudes{i} = [];
    IndividualBoutSyllEntropy{i} = [];
    IndividualBoutSyllMeanFreq{i} = [];
    BoutLens = BirdParameters(i).Bouts(:,2) - BirdParameters(i).Bouts(:,1);
    for j = 1:size(BirdParameters(i).Bouts, 1),
        MaxBoutLen = max(BoutLens);
        IndividualBoutSyllDurations{i}(end+1,:) = [BirdParameters(i).SAPFeatsMatrix(BirdParameters(i).Bouts(j,1):BirdParameters(i).Bouts(j,2),1)' (ones(1, MaxBoutLen - BoutLens(j)) * NaN)];
        IndividualBoutSyllAmplitudes{i}(end+1,:) = [BirdParameters(i).SAPFeatsMatrix(BirdParameters(i).Bouts(j,1):BirdParameters(i).Bouts(j,2),2)' (ones(1, MaxBoutLen - BoutLens(j)) * NaN)];
        IndividualBoutSyllEntropy{i}(end+1,:) = [BirdParameters(i).SAPFeatsMatrix(BirdParameters(i).Bouts(j,1):BirdParameters(i).Bouts(j,2),3)' (ones(1, MaxBoutLen - BoutLens(j)) * NaN)];
        IndividualBoutSyllMeanFreq{i}(end+1,:) = [BirdParameters(i).SAPFeatsMatrix(BirdParameters(i).Bouts(j,1):BirdParameters(i).Bouts(j,2),4)' (ones(1, MaxBoutLen - BoutLens(j)) * NaN)];
    end
    SyllableDistributionBins{i} = 0:0.005:max(IndividualBoutSyllDurations{i}(:));
    NotNaNIndices = ~isnan(IndividualBoutSyllDurations{i}(:));
    FullSyllDurHist{i}(1,:) = histc(IndividualBoutSyllDurations{i}(NotNaNIndices), SyllableDistributionBins{i})/length(NotNaNIndices);
    for j = 1:size(IndividualBoutSyllDurations{i},2),
        Not_NaN_Indices = ~isnan(IndividualBoutSyllDurations{i}(:,j));
        TempDurHist = histc(IndividualBoutSyllDurations{i}(Not_NaN_Indices,j), SyllableDistributionBins{i})/length(Not_NaN_Indices);
        if (~isempty(TempDurHist))
            SyllDurHist{i}(j,:) = TempDurHist(:)';
        else
            SyllDurHist{i}(j,:) = ones(size(SyllableDistributionBins{i})) * NaN;
        end
        MeanIndividualBoutSyllDurations{i}(j,:) = [j mean(IndividualBoutSyllDurations{i}(Not_NaN_Indices,j)) std(IndividualBoutSyllDurations{i}(Not_NaN_Indices,j)) length(find(Not_NaN_Indices))];
        MeanIndividualBoutSyllAmplitudes{i}(j,:) = [j mean(IndividualBoutSyllAmplitudes{i}(Not_NaN_Indices,j)) std(IndividualBoutSyllAmplitudes{i}(Not_NaN_Indices,j))];
        MeanIndividualBoutSyllEntropy{i}(j,:) = [j mean(IndividualBoutSyllEntropy{i}(Not_NaN_Indices,j)) std(IndividualBoutSyllEntropy{i}(Not_NaN_Indices,j))];
        MeanIndividualBoutSyllMeanFreq{i}(j,:) = [j mean(IndividualBoutSyllMeanFreq{i}(Not_NaN_Indices,j)) std(IndividualBoutSyllMeanFreq{i}(Not_NaN_Indices,j))];
        CorrIndividualBoutSyllDurHist{i}(j) = (SyllDurHist{i}(j,:) * FullSyllDurHist{i}')/(norm(FullSyllDurHist{i})*norm(SyllDurHist{i}(j,:)));
    end
    NotNaNIndices = ~isnan(IndividualBoutSyllDurations{i}(:));
    figure;
    % Plot # of syllables at each position
    plot(MeanIndividualBoutSyllDurations{i}(:,end), 'ko-');
    axis tight;
    Temp = axis;
    Temp = [0 MeanIndividualBoutSyllDurations{i}(end,1) 0 1.05*Temp(4)];
    axis(Temp);
    xlabel('Position in the bout', 'FontSize', 16);
    ylabel('# of syllables', 'FontSize', 16);
    
    % Plot mean and std for syllables at each position
    figure;
    errorbar(MeanIndividualBoutSyllDurations{i}(:,2), MeanIndividualBoutSyllDurations{i}(:,3), 'ko-');
    hold on;
    errorbar(0, mean(IndividualBoutSyllDurations{i}(NotNaNIndices)), std(IndividualBoutSyllDurations{i}(NotNaNIndices)), 'ks', 'MarkerSize', 6, 'LineWidth', 2);
    axis tight;
    Temp = axis;
    Temp = [-1 MeanIndividualBoutSyllDurations{i}(end,1) 0 1.05*Temp(4)];
    axis(Temp);
    xlabel('Position in the bout', 'FontSize', 16);
    ylabel('Mean syllable duration (sec)', 'FontSize', 16);
    
    % Plot examples of syllable duration histograms
    figure;
    for j = 1:10,
        subplot(5,2,j);
        plot(SyllableDistributionBins{i}, FullSyllDurHist{i}*100, 'r', 'LineWidth', 2);
        hold on;
        plot(SyllableDistributionBins{i}, SyllDurHist{i}(j,:)*100, 'b', 'LineWidth', 2);
        axis tight;
        Temp = axis;
        Temp = [0 SyllableDistributionBins{i}(end) 0 1.05*Temp(4)];
        axis(Temp);
        title(['Position #', num2str(j)], 'FontSize', 16);
        if (j > 8)
            xlabel('Time (sec)', 'FontSize', 16);
        end
        if (j == 5)
            ylabel('% of syllables', 'FontSize', 16);
        end
    end
    
    % Plot CV of syllable durations
    figure;
    plot(MeanIndividualBoutSyllDurations{1}(:,3)./MeanIndividualBoutSyllDurations{1}(:,2), 'ko-');
    hold on;
    NotNaNIndices = ~isnan(IndividualBoutSyllDurations{i}(:));
    plot([0 MeanIndividualBoutSyllDurations{1}(end,1)], ones(1,2)*(std(IndividualBoutSyllDurations{i}(find(NotNaNIndices))/mean(IndividualBoutSyllDurations{i}(find(NotNaNIndices))))), 'k--');
    axis tight;
    Temp = axis;
    Temp = [0 MeanIndividualBoutSyllDurations{1}(end,1) 0 1.05*Temp(4)];
    axis(Temp);
    xlabel('Position in the bout', 'FontSize', 16);
    ylabel('CV of syllable duration', 'FontSize', 16);
    
    % Plot mean and std for syllable amplitudes at each position
    figure;
    errorbar(MeanIndividualBoutSyllAmplitudes{i}(:,2), MeanIndividualBoutSyllAmplitudes{i}(:,3), 'ko-');
    hold on;
    errorbar(0, mean(IndividualBoutSyllAmplitudes{i}(NotNaNIndices)), std(IndividualBoutSyllAmplitudes{i}(NotNaNIndices)), 'ks', 'MarkerSize', 6, 'LineWidth', 2);
    axis tight;
    Temp = axis;
    Temp = [-1 MeanIndividualBoutSyllDurations{1}(end,1) 0 1.05*Temp(4)];
    axis(Temp);
    xlabel('Position in the bout', 'FontSize', 16);
    ylabel('Mean syllable amplitude (db)', 'FontSize', 16);
    
    % Plot CV of syllable amplitudes
    figure;
    plot(MeanIndividualBoutSyllAmplitudes{1}(:,3)./MeanIndividualBoutSyllAmplitudes{1}(:,2), 'ko-');
    hold on;
    NotNaNIndices = ~isnan(IndividualBoutSyllDurations{i}(:));
    plot([0 MeanIndividualBoutSyllDurations{1}(end,1)], ones(1,2)*(std(IndividualBoutSyllAmplitudes{i}(find(NotNaNIndices))/mean(IndividualBoutSyllAmplitudes{i}(find(NotNaNIndices))))), 'k--');
    axis tight;
    Temp = axis;
    Temp = [0 MeanIndividualBoutSyllDurations{1}(end,1) 0 1.05*Temp(4)];
    axis(Temp);
    xlabel('Position in the bout', 'FontSize', 16);
    ylabel('CV of syllable amplitude', 'FontSize', 16);
    
    % Plot mean and std for syllable entropy at each position
    figure;
    errorbar(MeanIndividualBoutSyllEntropy{i}(:,2), MeanIndividualBoutSyllEntropy{i}(:,3), 'ko-');
    hold on;
    errorbar(0, mean(IndividualBoutSyllEntropy{i}(NotNaNIndices)), std(IndividualBoutSyllEntropy{i}(NotNaNIndices)), 'ks', 'MarkerSize', 6, 'LineWidth', 2);
    axis tight;
    Temp = axis;
    Temp = [-1 MeanIndividualBoutSyllDurations{1}(end,1) 1.05*Temp(3) 0.95*Temp(4)];
    axis(Temp);
    xlabel('Position in the bout', 'FontSize', 16);
    ylabel('Mean syllable entropy', 'FontSize', 16);
    
    % Plot CV of syllable entropy
    figure;
    plot(abs(MeanIndividualBoutSyllEntropy{1}(:,3)./MeanIndividualBoutSyllEntropy{1}(:,2)), 'ko-');
    hold on;
    NotNaNIndices = ~isnan(IndividualBoutSyllDurations{i}(:));
    plot([0 MeanIndividualBoutSyllDurations{1}(end,1)], ones(1,2)*(std(IndividualBoutSyllEntropy{i}(find(NotNaNIndices))/mean(IndividualBoutSyllEntropy{i}(find(NotNaNIndices))))), 'k--');
    axis tight;
    Temp = axis;
    Temp = [0 MeanIndividualBoutSyllDurations{1}(end,1) 0 1.05*Temp(4)];
    axis(Temp);
    xlabel('Position in the bout', 'FontSize', 16);
    ylabel('CV of syllable entropy', 'FontSize', 16);
    
    % Plot mean and std for syllable mean frequencies at each position
    figure;
    errorbar(MeanIndividualBoutSyllMeanFreq{i}(:,2), MeanIndividualBoutSyllMeanFreq{i}(:,3), 'ko-');
    hold on;
    errorbar(0, mean(IndividualBoutSyllMeanFreq{i}(NotNaNIndices)), std(IndividualBoutSyllMeanFreq{i}(NotNaNIndices)), 'ks', 'MarkerSize', 6, 'LineWidth', 2);
    axis tight;
    Temp = axis;
    Temp = [-1 MeanIndividualBoutSyllDurations{1}(end,1) 0 1.05*Temp(4)];
    axis(Temp);
    xlabel('Position in the bout', 'FontSize', 16);
    ylabel('Mean syllable frequency (Hz)', 'FontSize', 16);
    
    % Plot CV of syllable mean frequency
    figure;
    plot(MeanIndividualBoutSyllMeanFreq{1}(:,3)./MeanIndividualBoutSyllMeanFreq{1}(:,2), 'ko-');
    hold on;
    NotNaNIndices = ~isnan(IndividualBoutSyllDurations{i}(:));
    plot([0 MeanIndividualBoutSyllDurations{1}(end,1)], ones(1,2)*(std(IndividualBoutSyllMeanFreq{i}(find(NotNaNIndices))/mean(IndividualBoutSyllMeanFreq{i}(find(NotNaNIndices))))), 'k--');
    axis tight;
    Temp = axis;
    Temp = [0 MeanIndividualBoutSyllDurations{1}(end,1) 0 1.05*Temp(4)];
    axis(Temp);
    xlabel('Position in the bout', 'FontSize', 16);
    ylabel('CV of syllable frequency', 'FontSize', 16);
    
end