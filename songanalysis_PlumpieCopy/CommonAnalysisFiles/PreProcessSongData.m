function [] = PreProcessSongData(BirdDetailsTextFile, MinTrialNo)

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
        MeanIndividualBoutSyllDurations{i}(j,:) = [j mean(IndividualBoutSyllDurations{i}(Not_NaN_Indices,j)) std(IndividualBoutSyllDurations{i}(Not_NaN_Indices,j))];
        MeanIndividualBoutSyllAmplitudes{i}(j,:) = [j mean(IndividualBoutSyllAmplitudes{i}(Not_NaN_Indices,j)) std(IndividualBoutSyllAmplitudes{i}(Not_NaN_Indices,j))];
        MeanIndividualBoutSyllEntropy{i}(j,:) = [j mean(IndividualBoutSyllEntropy{i}(Not_NaN_Indices,j)) std(IndividualBoutSyllEntropy{i}(Not_NaN_Indices,j))];
        MeanIndividualBoutSyllMeanFreq{i}(j,:) = [j mean(IndividualBoutSyllMeanFreq{i}(Not_NaN_Indices,j)) std(IndividualBoutSyllMeanFreq{i}(Not_NaN_Indices,j))];
        CorrIndividualBoutSyllDurHist{i}(j) = (SyllDurHist{i}(j,:) * FullSyllDurHist{i}')/(norm(FullSyllDurHist{i})*norm(SyllDurHist{i}(j,:)));
    end
end
% Now to find the INs in this
% First I will just label all motifs


% Now collate all labels into bouts
disp('Collating all labels ...');
for i = 1:length(BirdParameters),
    [BirdParameters(i).BoutLabels] = LSINA_CollateBouts(BirdParameters(i));
    if (BirdParameters(i).Continuousdata == 1)
        BirdParameters(i).ActualAllLabels = BirdParameters(i).AllLabels;
        BirdParameters(i).ActualAllOnsets = BirdParameters(i).AllOnsets;
        BirdParameters(i).ActualAllOffsets = BirdParameters(i).AllOffsets;
        
        Indices = find((BirdParameters(i).BoutLabels ~= 'Q') & (BirdParameters(i).BoutLabels ~= 'q'));
        
        BirdParameters(i).AllLabels = BirdParameters(i).BoutLabels(Indices);
        BirdParameters(i).AllOnsets = [];
        BirdParameters(i).AllOffsets = [];
        for j = 1:size(BirdParameters(i).BoutIndices,1),
            BirdParameters(i).AllOnsets = [BirdParameters(i).AllOnsets; (BirdParameters(i).ActualAllOnsets(BirdParameters(i).BoutIndices(j,1):BirdParameters(i).BoutIndices(j,2)))];
            BirdParameters(i).AllOffsets = [BirdParameters(i).AllOffsets; (BirdParameters(i).ActualAllOffsets(BirdParameters(i).BoutIndices(j,1):BirdParameters(i).BoutIndices(j,2)))];
        end
    end
end

% Now calculate SAP features for all the files
disp('Calculating SAP features ...');
for i = 1:length(BirdParameters),
    disp(['Calculating SAP features for ', BirdParameters(i).BirdName, '-', BirdParameters(i).DataLabel, ' ...']);
    if (exist([BirdParameters(i).BirdName, '.', BirdParameters(i).DataLabel, '.SAPFeats.mat'], 'file'))
        load([BirdParameters(i).BirdName, '.', BirdParameters(i).DataLabel, '.SAPFeats.mat']);
        BirdParameters(i).SAPFeats = SAPFeats;
    else
        [BirdParameters(i).SAPFeats] = LSINA_CalcSAPFeats(BirdParameters(i));
        SAPFeats = BirdParameters(i).SAPFeats;
        save([BirdParameters(i).BirdName, '.', BirdParameters(i).DataLabel, '.SAPFeats.mat'], 'SAPFeats');
    end
end

% Now find out all the bout starting syllables, the first occuring motif
% syllable in a bout and the syllable just before that
disp('Characterizing bout starting ...');

[UniqueBirds, UniqueBirdIndices] = unique([BirdParameters.SerialNo]);

for i = 1:length(BirdParameters),
    [BirdParameters(i).BoutStartSyll, BirdParameters(i).FirstMotifSyll, BirdParameters(i).PreFirstMotifSyll, BirdParameters(i).FirstMotifSyllIndex] = LSINA_CalculateBoutStartingSyllStats(BirdParameters(i));
    
    BirdParameters(i).DiffBoutStartSylls = unique(BirdParameters(i).BoutStartSyll);
    BirdParameters(i).NumDiffBoutStartSylls = length(BirdParameters(i).DiffBoutStartSylls);
    for j = 1:BirdParameters(i).NumDiffBoutStartSylls,
        BirdParameters(i).ProportionBoutStartSylls(j) = length(find(BirdParameters(i).BoutStartSyll == BirdParameters(i).DiffBoutStartSylls(j)))/length(BirdParameters(i).BoutStartSyll);
    end
    
    BirdParameters(i).DiffFirstMotifSylls = unique(BirdParameters(i).FirstMotifSyll);    
    BirdParameters(i).NumDiffFirstMotifSylls = length(BirdParameters(i).DiffFirstMotifSylls);
    for j = 1:BirdParameters(i).NumDiffFirstMotifSylls,
        BirdParameters(i).ProportionFirstMotifSylls(j) = length(find(BirdParameters(i).FirstMotifSyll == BirdParameters(i).DiffFirstMotifSylls(j)))/length(BirdParameters(i).FirstMotifSyll);
    end
    
    BirdParameters(i).DiffPreFirstMotifSylls = unique(BirdParameters(i).PreFirstMotifSyll);
    BirdParameters(i).NumDiffPreFirstMotifSylls = length(BirdParameters(i).DiffPreFirstMotifSylls);
    for j = 1:BirdParameters(i).NumDiffPreFirstMotifSylls,
        BirdParameters(i).ProportionPreFirstMotifSylls(j) = length(find(BirdParameters(i).PreFirstMotifSyll == BirdParameters(i).DiffPreFirstMotifSylls(j)))/length(BirdParameters(i).PreFirstMotifSyll);
    end
end

% Now identify INs and put position labels for all of the INs
disp('IN identification ...');
for i = 1:length(BirdParameters),
    [BirdParameters(i).INResults] = LSINA_IdentifyINs(BirdParameters(i));
end
    
% Write the data to a text file - the labels, the motif postions, syllable
% positions and IN positions for verification
for i = 1:length(BirdParameters),
    Fid = fopen([BirdParameters(i).BirdName, '.', BirdParameters(i).DataLabel, '.Labels.txt'], 'w');
    fprintf(Fid, 'Actual Labels\tMotif Position\tSyllable Position\tIN Position\n');
    for j = 1:length(BirdParameters(i).AllLabels),
        fprintf(Fid, '%c\t%d\t%d\t%d\t%d\n', BirdParameters(i).AllLabels(j), BirdParameters(i).INResults.AllBoutStartINs(2,j), BirdParameters(i).INResults.AllCommonMotifPositions(1,j), BirdParameters(i).INResults.AllCommonMotifPreSyllables(end,j),BirdParameters(i).INResults.AllCommonMotifPreINs(end,j));
    end
    fclose(Fid);
end

% Now one more thing to do, I need to have a vector that tells me which
% bout comes from which file. Bout Indices is useful for this in
% non-continuous data

for i = 1:length(BirdParameters),
    if (BirdParameters(i).Continuousdata == 0)
        Bout_Files = cumsum(cellfun(@size, BirdParameters(i).BoutIndices, num2cell(ones(size(BirdParameters(i).BoutIndices)))));
        for j = 1:length(Bout_Files),
            if (j == 1)
                BirdParameters(i).Bout_Files(1:Bout_Files(j)) = j;
            else
                BirdParameters(i).Bout_Files((Bout_Files(j-1)+1):Bout_Files(j)) = j;
            end
        end
    else
        FileTimes = cumsum(BirdParameters(i).FileLen);
        BoutOnsetTimes = BirdParameters(i).ActualAllOnsets(BirdParameters(i).BoutIndices(:,1));
        BoutOffsetTimes = BirdParameters(i).ActualAllOffsets(BirdParameters(i).BoutIndices(:,2));
        for k = 1:length(BoutOnsetTimes),
            BirdParameters(i).Bout_Files(k,1) = find(FileTimes >= BoutOnsetTimes(k), 1, 'first');
            BirdParameters(i).Bout_Files(k,2) = find(FileTimes >= BoutOffsetTimes(k), 1, 'first');
        end
    end
end

% Now put all the data in one common format called Syllable Data. This will
% be another field within BirdParameters - this will contain all the data
% and will correspond to the labels in BirdParameters.AllLabels

for i = 1:length(BirdParameters),
%    BirdParameters(i).SyllableData(:,1) =  
end
disp('Finished pre-processing data');