function [] = AnalyzeContinuousData(BirdDetailsTextFile, MinTrialNo)

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
    [BirdParameters(i).SongFileNames] =LSINA_GetDataFileNames(BirdParameters(i));
end

% Now load up the note files and the length of each file
disp('Loading up note data ...');
for i = 1:length(BirdParameters),
    [BirdParameters(i).NoteInfo, BirdParameters(i).FileLen] = LSINA_LoadNoteFileInfo(BirdParameters(i));
end

% Now split up the files into bouts based on inter-bout interval that is
% also specified in the .csv file
disp('Identifying bouts ...');
for i = 1:length(BirdParameters),
    [BirdParameters(i).Bouts, BirdParameters(i).BoutIndices, BirdParameters(i).AllOnsets, BirdParameters(i).AllOffsets, BirdParameters(i).AllLabels] = LSINA_IdentifyBouts(BirdParameters(i));
end

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

disp('Finished');