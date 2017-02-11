function [] = PreProcessSongDataIntoBouts(BirdDetailsTextFile)

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
    [BirdParameters(i).SongFileNames] = LSINA_GetDataFileNames(BirdParameters(i));
end

% Now load up the note files and the length of each file
disp('Loading up note data ...');
for i = 1:length(BirdParameters),
    [BirdParameters(i).NoteInfo, BirdParameters(i).FileLen] = LSINA_LoadNoteFileInfo(BirdParameters(i));
end

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

% =========================================================================
% Plots
% Bout length for different birds - mean and variability of bout length
for i = 1:length(BirdParameters),
    if (BirdParameters(i).Continuousdata == 0)
        for j = 1:length(BirdParameters(i).BoutIndices),
            for k = 1:size(BirdParameters(i).BoutIndices{j},1),
                BirdParameters(i).BoutLen{j}(k) = BirdParameters(i).NoteInfo{j}.offsets(BirdParameters(i).BoutIndices{j}(k,2)) - BirdParameters(i).NoteInfo{j}.onsets(BirdParameters(i).BoutIndices{j}(k,1));
            end
        end
        BirdParameters(i).MeanBoutLen = mean(cell2mat(BirdParameters(i).BoutLen));
        BirdParameters(i).STDBoutLen = std(cell2mat(BirdParameters(i).BoutLen));
        BirdParameters(i).SEMBoutLen = BirdParameters(i).STDBoutLen/sqrt(length(cell2mat(BirdParameters(i).BoutLen)));
        BirdParameters(i).NumBouts = length(cell2mat(BirdParameters(i).BoutLen));
    else
        for j = 1:size(BirdParameters(i).BoutIndices,1),
            BirdParameters(i).BoutLen(j) = BirdParameters(i).ActualAllOffsets(BirdParameters(i).BoutIndices(j,2)) - BirdParameters(i).ActualAllOnsets(BirdParameters(i).BoutIndices(j,1));
        end
        BirdParameters(i).MeanBoutLen = mean(BirdParameters(i).BoutLen);
        BirdParameters(i).STDBoutLen = std(BirdParameters(i).BoutLen);
        BirdParameters(i).SEMBoutLen = BirdParameters(i).STDBoutLen/sqrt(length(BirdParameters(i).BoutLen));        
        BirdParameters(i).NumBouts = length(BirdParameters(i).BoutLen);
    end
    XAxisLabelStrings{i} = [BirdParameters(i).BirdName, '.', BirdParameters(i).DataLabel];
end

% For birds that have multiple days I'm choosing the day with maximum
% number of bouts
[UniqueBirds, UniqueBirdIndices] = unique([BirdParameters.SerialNo]);
for i = 1:length(UniqueBirdIndices),
    BirdDays = find([BirdParameters.SerialNo] == i);
    [MaxBouts, DayWithMaxBouts_Index] = max([BirdParameters(BirdDays).NumBouts]);
    UniqueBirdIndices(i) = BirdDays(DayWithMaxBouts_Index);
end

% First plot mean of bout length for all birds
figure;
MeanBoutLenPlots = errorbar([BirdParameters.MeanBoutLen], [BirdParameters.STDBoutLen], 'ks');
axis tight;
PlotAxis = axis;
PlotAxis = [0.5 (length(BirdParameters) + 0.5) 0 1.05*PlotAxis(4)];
for i = 1:length(BirdParameters),
    text(i-0.4, 0.85*(BirdParameters(i).MeanBoutLen - BirdParameters(i).STDBoutLen), ['(', num2str(BirdParameters(i).NumBouts), ')']);
end
LSINA_FixPlotLabels(gcf, 8, 'Fig_Position', [94 195 1100 350], 'Fig_Axis', PlotAxis, 'Fig_XTickLabels', XAxisLabelStrings, 'Fig_XTickLabel_Rotation', 45, 'Fig_Title', 'Mean +/- STD Bout Length across birds', 'XAxis_Label', 'Bird #', 'YAxis_Label', 'Bout length (ms)', 'Fig_Save', 'MeanBoutLen.png');

% Next plot CV of bout length for all birds
figure;
CVBoutLenPlots = plot([BirdParameters.STDBoutLen]./[BirdParameters.MeanBoutLen], 'ks');
axis tight;
PlotAxis = axis;
PlotAxis = [0.5 (length(BirdParameters) + 0.5) 0 1.05*PlotAxis(4)];
for i = 1:length(BirdParameters),
    text(i-0.4, 0.75*(BirdParameters(i).STDBoutLen/BirdParameters(i).MeanBoutLen), ['(', num2str(BirdParameters(i).NumBouts), ')']);
end
LSINA_FixPlotLabels(gcf, 10, 'Fig_Position', [94 195 1100 350], 'Fig_Axis', PlotAxis, 'Fig_XTickLabels', XAxisLabelStrings, 'Fig_XTickLabel_Rotation', 45, 'Fig_Title', 'CV (std/mean) of Bout Length across birds', 'XAxis_Label', 'Bird #', 'YAxis_Label', 'CV', 'Fig_Save', 'CVBoutLen.png');

% Plots of number of motifs per bout


% 1. The # of different syllables that are seen at the start of a bout of
% song for each of the birds
% 2. The # of different syllables that occur before the first motif
% syllable
% 3. The # of different first motif syllable in a bout of song

% Individual bird plots
figure;
BoutBeginningPlots = bar([[BirdParameters(UniqueBirdIndices).NumDiffBoutStartSylls]' [BirdParameters(UniqueBirdIndices).NumDiffPreFirstMotifSylls]' [BirdParameters(UniqueBirdIndices).NumDiffFirstMotifSylls]']);
LSINA_FixPlotLabels(gcf, 14, 'Fig_Title', 'Characterization of bout initiation syllables - individual birds', 'XAxis_Label', 'Bird #', 'YAxis_Label', '#', 'Fig_Legend', [{'Bout Start Syllables'} {'Pre-Motif start syllables'} {'First motif syllables'}], 'Fig_Save', 'BoutStartingCharacterization.png');

% Summary across all birds
figure;
errorbar(1, mean([BirdParameters(UniqueBirdIndices).NumDiffBoutStartSylls]), std([BirdParameters(UniqueBirdIndices).NumDiffBoutStartSylls])/length(UniqueBirdIndices), 'ks', 'LineWidth', 1, 'MarkerSize', 10);
hold on;
errorbar(2, mean([BirdParameters(UniqueBirdIndices).NumDiffPreFirstMotifSylls]), std([BirdParameters(UniqueBirdIndices).NumDiffPreFirstMotifSylls])/length(UniqueBirdIndices), 'ks', 'LineWidth', 1, 'MarkerSize', 10);
errorbar(3, mean([BirdParameters(UniqueBirdIndices).NumDiffFirstMotifSylls]), std([BirdParameters(UniqueBirdIndices).NumDiffFirstMotifSylls])/length(UniqueBirdIndices), 'ks', 'LineWidth', 1, 'MarkerSize', 10);
axis tight;
Temp = axis;

PlotAxis = [0.5 3.5 0 1.1*Temp(4)];

LSINA_FixPlotLabels(gcf, 14, 'Fig_Title', ['Characterization of bout initiation syllables - summary across birds (n=', num2str(length(UniqueBirdIndices)), ')'], 'YAxis_Label', '#', 'Fig_XTickLabels', [{'Bout Start Syllables'} {'Pre-Motif start syllables'} {'First motif syllables'}], 'Fig_Axis', PlotAxis, 'Fig_Save', 'BoutStartingCharacterization_Summary.png');

clear Temp;
% Now plot the diversity and proportion of bout initiation syllables
for i = 1:length(BirdParameters),
    Temp(i).NumDiffSylls = BirdParameters(i).NumDiffBoutStartSylls;
    Temp(i).DiffSylls = BirdParameters(i).DiffBoutStartSylls;
    Temp(i).ProportionSylls = BirdParameters(i).ProportionBoutStartSylls;
    Temp(i).BirdName = BirdParameters(i).BirdName;
    Temp(i).DataLabel = BirdParameters(i).DataLabel;
end

LSINA_PlotProportionDiversitySyllables(Temp, 'Bout Start', 'BoutStartingDiversityProportions.png');

clear Temp;
% Now plot the diversity and proportion of first motif syllables
for i = 1:length(BirdParameters),
    Temp(i).NumDiffSylls = BirdParameters(i).NumDiffFirstMotifSylls;
    Temp(i).DiffSylls = BirdParameters(i).DiffFirstMotifSylls;
    Temp(i).ProportionSylls = BirdParameters(i).ProportionFirstMotifSylls;
    Temp(i).BirdName = BirdParameters(i).BirdName;
    Temp(i).DataLabel = BirdParameters(i).DataLabel;
end

LSINA_PlotProportionDiversitySyllables(Temp, 'First Motif', 'FirstMotifDiversityProportions.png');

clear Temp;
% Now plot the diversity and proportion of pre first motif syllables
for i = 1:length(BirdParameters),
    Temp(i).NumDiffSylls = BirdParameters(i).NumDiffPreFirstMotifSylls;
    Temp(i).DiffSylls = BirdParameters(i).DiffPreFirstMotifSylls;
    Temp(i).ProportionSylls = BirdParameters(i).ProportionPreFirstMotifSylls;
    Temp(i).BirdName = BirdParameters(i).BirdName;
    Temp(i).DataLabel = BirdParameters(i).DataLabel;
end

LSINA_PlotProportionDiversitySyllables(Temp, 'Pre First Motif', 'PreFirstMotifDiversityProportions.png');

% Now plot # of INs for each bird if it starts the bout with only INs and
% if it starts with something else
% I'm going to do this individually for each bird first and then I can
% generate a summary plot
Colors = 'rbcmyk';

for i = 1:length(UniqueBirdIndices),
    FigLegend = [];
    figure;
    Indices = [];
    DiffBoutStartSylls = [];
    NumTrials = [];
    for j = 1:length(BirdParameters),
        if (~isempty(strfind(BirdParameters(j).BirdName, BirdParameters(UniqueBirdIndices(i)).BirdName)))
            Indices(end+1) = j;
            FigLegend{end+1} = BirdParameters(j).DataLabel;
            DiffBoutStartSylls = [DiffBoutStartSylls; BirdParameters(j).DiffBoutStartSylls(:)];
            NumTrials = [NumTrials; (BirdParameters(j).ProportionBoutStartSylls(:) * length(BirdParameters(j).INResults.NumBoutStartINs))];
        end
    end
    [DiffBoutStartSylls, UniqueIndices] = unique(DiffBoutStartSylls);
    NumTrials = NumTrials(UniqueIndices);
    
    DiffBoutStartSylls = DiffBoutStartSylls(find(NumTrials >= MinTrialNo));
    MeanNumINs = [];
    StdNumINs = [];
    SEMNumINs = [];
    for j = 1:length(DiffBoutStartSylls),
        for k = 1:length(Indices),
            MeanNumINs(k,j) = mean(BirdParameters(Indices(k)).INResults.NumBoutStartINs(find(BirdParameters(Indices(k)).BoutStartSyll == DiffBoutStartSylls(j))));
            STDNumINs(k,j) = std(BirdParameters(Indices(k)).INResults.NumBoutStartINs(find(BirdParameters(Indices(k)).BoutStartSyll == DiffBoutStartSylls(j))));
            SEMNumINs(k,j) = STDNumINs(k,j)/sqrt(length(find(BirdParameters(Indices(k)).BoutStartSyll == DiffBoutStartSylls(j))));
            TrialNum(k,j) = length(find(BirdParameters(Indices(k)).BoutStartSyll == DiffBoutStartSylls(j)));
        end
    end
    for k = 1:size(MeanNumINs,1),
        errorbar(MeanNumINs(k,:), SEMNumINs(k,:), [Colors(k),'o-']);
        hold on;
        for j = 1:size(MeanNumINs, 2),
            text(j, (0.5 + (k-1)*0.2), ['(', num2str(TrialNum(k,j)), ')'], 'Color', Colors(k));
        end
    end
    axis tight;
    Temp = axis;
    PlotAxis = [0.75 (length(DiffBoutStartSylls) + 1.5) 0 1.1*Temp(4)];
    XAxis_LabelStrings = mat2cell(DiffBoutStartSylls, ones(size(DiffBoutStartSylls,1),1), ones(size(DiffBoutStartSylls,2),1));
    LSINA_FixPlotLabels(gcf, 14, 'Fig_Position', [94 195 450 300], 'Fig_Legend', FigLegend, 'Fig_Title', [BirdParameters(UniqueBirdIndices(i)).BirdName, ': # of INs: different starts'], 'YAxis_Label', 'Mean # of INs', 'Fig_XTickLabels', XAxis_LabelStrings, 'Fig_Axis', PlotAxis, 'Fig_Save', [BirdParameters(UniqueBirdIndices(i)).BirdName, '.NumINs.png']);
end

% Now plot the mean # of INs for each bird individually depending on
% whether the syllables at the start of the bout were only INs, or they
% contained a diverse combination of other syllables as well.

for i = 1:length(UniqueBirdIndices),
    FigLegend = [];
    figure;
    Indices = [];
    DiffBoutStartSylls = [];
    NumTrials = [];
    for j = 1:length(BirdParameters),
        if (~isempty(strfind(BirdParameters(j).BirdName, BirdParameters(UniqueBirdIndices(i)).BirdName)))
            Indices(end+1) = j;
            FigLegend{end+1} = BirdParameters(j).DataLabel;
        end
    end
    
    for j = 1:length(Indices),
        TotalSyllCount = 0;
        NumINs_INBouts = [];
        NumINs_NonINBouts = [];
        for k = 1:length(BirdParameters(Indices(j)).INResults.IndividualBoutLabels),
            FirstMotifSyll = find((BirdParameters(Indices(j)).INResults.AllCommonMotifPositions(1,:) > 0) & (BirdParameters(Indices(j)).INResults.AllCommonMotifPositions(2,:) == k), 1, 'first');
            FirstMotifSyll = FirstMotifSyll - TotalSyllCount;
            if ((FirstMotifSyll - 1) == sum(BirdParameters(Indices(j)).INResults.INs{k}(1:FirstMotifSyll-1)))
                NumINs_INBouts(end+1) = BirdParameters(Indices(j)).INResults.NumBoutStartINs(k);
                BirdParameters(Indices(j)).IN_NonIN_BoutIndices(k) = 1;
            else
                NumINs_NonINBouts(end+1) = BirdParameters(Indices(j)).INResults.NumBoutStartINs(k);
                BirdParameters(Indices(j)).IN_NonIN_BoutIndices(k) = 0;
            end
            TotalSyllCount = TotalSyllCount + length(BirdParameters(Indices(j)).INResults.IndividualBoutLabels{k});
        end
        
        errorbar([mean(NumINs_INBouts) mean(NumINs_NonINBouts)], [std(NumINs_INBouts)/sqrt(length(NumINs_INBouts)) std(NumINs_NonINBouts)/sqrt(length(NumINs_NonINBouts))], [Colors(j),'o-']);
        hold on;
        text(1, (0.5 + (j-1)*0.2), ['(', num2str(length(NumINs_INBouts)), ')'], 'Color', Colors(j));
        text(2, (0.5 + (j-1)*0.2), ['(', num2str(length(NumINs_NonINBouts)), ')'], 'Color', Colors(j));
    end
    axis tight;
    Temp = axis;
    PlotAxis = [0.5 3 0 1.1*Temp(4)];
    XAxis_LabelStrings = [{'IN only bouts'} {'Other bouts'}];
    LSINA_FixPlotLabels(gcf, 14, 'Fig_Position', [94 195 700 300], 'Fig_Title', [BirdParameters(UniqueBirdIndices(i)).BirdName, ': # of INs: IN only starts and other starts'], 'YAxis_Label', 'Mean # of INs', 'Fig_XTickLabels', XAxis_LabelStrings, 'Fig_Legend', FigLegend, 'Fig_Axis', PlotAxis, 'Fig_Save', [BirdParameters(UniqueBirdIndices(i)).BirdName, '.NumINs.INBouts_NonINBouts.png']);
end

% Now plot the summary for all birds for the mean # of INs in bouts with
% only INs and bouts with other syllables too. For individual birds with
% multiple days, I'm only going to consider the first day of data

% First plot the proportion of bouts on any given day that have only INs
% and that have other syllables as well
for i = 1:length(BirdParameters),
    Fraction_INBouts(i) = length(find(BirdParameters(i).IN_NonIN_BoutIndices))/length(BirdParameters(i).IN_NonIN_BoutIndices);
end

PlotAxis = axis;
PlotAxis = [0.5 (length(BirdParameters) + 0.5) 0 1.05*PlotAxis(4)];
for i = 1:length(BirdParameters),
    text(i-0.4, 0.75*(BirdParameters(i).STDBoutLen/BirdParameters(i).MeanBoutLen), ['(', num2str(BirdParameters(i).NumBouts), ')']);
end
LSINA_FixPlotLabels(gcf, 10, 'Fig_Position', [94 195 1100 350], 'Fig_Axis', PlotAxis, 'Fig_XTickLabels', XAxisLabelStrings, 'Fig_XTickLabel_Rotation', 45, 'Fig_Title', 'CV (std/mean) of Bout Length across birds', 'XAxis_Label', 'Bird #', 'YAxis_Label', 'CV', 'Fig_Save', 'CVBoutLen.png');

clear MeanNumINs STDNumINs SEMNumINs;

for i = 1:length(UniqueBirdIndices),
    for j = 1:length(BirdParameters),
        if (~isempty(strfind(BirdParameters(j).BirdName, BirdParameters(UniqueBirdIndices(i)).BirdName)))
            Index = j;
            break;
        end
    end
    
    for j = Index,
        TotalSyllCount = 0;
        NumINs_INBouts = [];
        NumSylls_INBouts = [];
        TotalNumSylls_INBouts = [];
        NumINs_NonINBouts = [];
        NumSylls_NonINBouts = [];
        TotalNumSylls_NonINBouts = [];
        
        for k = 1:length(BirdParameters(j).INResults.IndividualBoutLabels),
            FirstMotifSyll = find((BirdParameters(j).INResults.AllCommonMotifPositions(1,:) > 0) & (BirdParameters(j).INResults.AllCommonMotifPositions(2,:) == k), 1, 'first');
            FirstMotifSyll = FirstMotifSyll - TotalSyllCount;
            if ((FirstMotifSyll - 1) == sum(BirdParameters(j).INResults.INs{k}(1:FirstMotifSyll-1)))
                NumINs_INBouts(end+1) = BirdParameters(j).INResults.NumBoutStartINs(k);
                NumSylls_INBouts(end+1) = BirdParameters(j).INResults.NumBoutStartSyllables(k);
                TotalNumSylls_INBouts(end+1) = BirdParameters(j).INResults.TotalNumBoutStartSyllables(k);
            else
                NumINs_NonINBouts(end+1) = BirdParameters(j).INResults.NumBoutStartINs(k);
                NumSylls_NonINBouts(end+1) = BirdParameters(j).INResults.NumBoutStartSyllables(k);
                TotalNumSylls_NonINBouts(end+1) = BirdParameters(j).INResults.TotalNumBoutStartSyllables(k);
            end
            TotalSyllCount = TotalSyllCount + length(BirdParameters(j).INResults.IndividualBoutLabels{k});
        end
        
        MeanNumINs(i,:) = [mean(NumINs_INBouts) mean(NumINs_NonINBouts)];
        STDNumINs(i,:) = [std(NumINs_INBouts) std(NumINs_NonINBouts)];
        SEMNumINs(i,:) = [(std(NumINs_INBouts)/sqrt(length(NumINs_INBouts))) (std(NumINs_NonINBouts)/sqrt(length(NumINs_NonINBouts)))];
        
        MeanNumSylls(i,:) = [mean(NumSylls_INBouts) mean(NumSylls_NonINBouts)];
        STDNumSylls(i,:) = [std(NumSylls_INBouts) std(NumSylls_NonINBouts)];
        SEMNumSylls(i,:) = [(std(NumSylls_INBouts)/sqrt(length(NumSylls_INBouts))) (std(NumSylls_NonINBouts)/sqrt(length(NumSylls_NonINBouts)))];
        
        MeanTotalNumSylls(i,:) = [mean(TotalNumSylls_INBouts) mean(TotalNumSylls_NonINBouts)];
        STDTotalNumSylls(i,:) = [std(TotalNumSylls_INBouts) std(TotalNumSylls_NonINBouts)];
        SEMTotalNumSylls(i,:) = [(std(TotalNumSylls_INBouts)/sqrt(length(TotalNumSylls_INBouts))) (std(TotalNumSylls_NonINBouts)/sqrt(length(TotalNumSylls_NonINBouts)))];
    end
end

figure;
hold on;
MeanINBouts_Bar = bar(mean(MeanNumINs));
set(MeanINBouts_Bar, 'FaceColor', 'none', 'EdgeColor', 'k');
for i = 1:size(MeanNumINs,2),
    errorbar(i, mean(MeanNumINs(:,i)), std(MeanNumINs(:,i))/sqrt(length(UniqueBirdIndices)), 'ks');
end
plot(repmat([1.1 1.9], length(UniqueBirdIndices), 1)', MeanNumINs', 'ko-');
plot([1.1 1.9], MeanNumINs(2,:), 'ro-');
axis tight;
Temp = axis;
PlotAxis = [0.5 2.5 0 1.1*Temp(4)];
XAxis_LabelStrings = [{'IN only bouts'} {'Other bouts'}];
LSINA_FixPlotLabels(gcf, 14, 'Fig_Position', [94 195 500 500], 'Fig_Title', '# of INs: IN only starts and other starts', 'YAxis_Label', 'Mean # of INs', 'Fig_XTickLabels', XAxis_LabelStrings, 'Fig_Axis', PlotAxis, 'Fig_Save', ['AllBirds.NumINs.INBouts_NonINBouts.png']);
p = signrank(MeanNumINs(:,1), MeanNumINs(:,2));
disp(['No of INs at the beginning of bouts with only INs vs. No of INs at the beginning of all other bouts: Wilcoxon sign-rank test p is ', num2str(p)]);


% Now plot the summary for all birds for the mean # of INs in bouts with
% only INs on successive days to assess short-term stability of # of INs.
% I'm going to consider anything < 1 week as short-term

clear MeanNumINs STDNumINs SEMNumINs;

for i = 1:length(UniqueBirdIndices),
    Indices = [];
    for j = 1:length(BirdParameters),
        if (~isempty(strfind(BirdParameters(j).BirdName, BirdParameters(UniqueBirdIndices(i)).BirdName)))
            Indices(end+1) = j;
            FigLegend{end+1} = BirdParameters(j).DataLabel;
        end
    end
    
    if (length(Indices) > 0)
        for j = 1:length(Indices),
            for k = 1:length(Indices),
                DaysBetweenRecordings{i}(j,k) = etime(datevec(BirdParameters(Indices(j)).Datestring, 'mmmm dd, yyyy'), datevec(BirdParameters(Indices(k)).Datestring, 'mmmm dd, yyyy'))/(24*3600);
                % etime gives the answer in seconds so I have divided by
                % (24*3600) to convert to days
            end
        end
    end
    
    for j = 1:length(Indices),
        TotalSyllCount = 0;
        NumINs_INBouts = [];
        for k = 1:length(BirdParameters(Indices(j)).INResults.IndividualBoutLabels),
            FirstMotifSyll = find((BirdParameters(Indices(j)).INResults.AllCommonMotifPositions(1,:) > 0) & (BirdParameters(Indices(j)).INResults.AllCommonMotifPositions(2,:) == k), 1, 'first');
            FirstMotifSyll = FirstMotifSyll - TotalSyllCount;
            if ((FirstMotifSyll - 1) == sum(BirdParameters(Indices(j)).INResults.INs{k}(1:FirstMotifSyll-1)))
                NumINs_INBouts(end+1) = BirdParameters(Indices(j)).INResults.NumBoutStartINs(k);
            end
            TotalSyllCount = TotalSyllCount + length(BirdParameters(Indices(j)).INResults.IndividualBoutLabels{k});
        end
        
        MeanNumINs{i}(j) = mean(NumINs_INBouts);
        STDNumINs{i}(j) = std(NumINs_INBouts);
        SEMNumINs{i}(j) = std(NumINs_INBouts)/sqrt(length(NumINs_INBouts));
    end
end

ShortTermMeanAcrossBirds = [];
LongTermMeanAcrossBirds = [];
for i = 1:length(DaysBetweenRecordings),
   [Day1, Day2] = find((DaysBetweenRecordings{i} < 3) & (DaysBetweenRecordings{i} > 0));
   if (~isempty(Day1))
       ShortTermMeanAcrossBirds(end+1,:) = [MeanNumINs{i}(Day2(1)) MeanNumINs{i}(Day1(1))];
   end
   
   [Day1, Day2] = find((DaysBetweenRecordings{i} >= 180));
   if (~isempty(Day1))
       LongTermMeanAcrossBirds(end+1,:) = [MeanNumINs{i}(Day2(1)) MeanNumINs{i}(Day1(1))];
   end
end

% Plot all birds with longitudinal data - difference in days on x-axis and
% mean # on y-axis
figure;
hold on;
for i = 1:length(DaysBetweenRecordings),
    if (length(DaysBetweenRecordings{i}) > 1)
        plot(DaysBetweenRecordings{i}(:,1), MeanNumINs{i}, 'ko-');
    end
end
axis tight;
Temp = axis;
PlotAxis = [0.5 1.1*Temp(2) 0 1.1*Temp(4)];
LSINA_FixPlotLabels(gcf, 14, 'Fig_Position', [94 195 500 500], 'Fig_Title', '# of INs on consecutive days', 'XAxis_Label', 'Days since 1st recording', 'YAxis_Label', 'Mean # of INs', 'Fig_Axis', PlotAxis, 'Fig_Save', ['AllBirds.LongitudinalData.NumINs.INBouts.png']);

% Plot summary across birds - short term
figure;
hold on;
ShortTermMeanAcrossBirds_Bar = bar(mean(ShortTermMeanAcrossBirds));
set(ShortTermMeanAcrossBirds_Bar, 'FaceColor', 'none', 'EdgeColor', 'k');
for i = 1:size(ShortTermMeanAcrossBirds,2),
    errorbar(i, mean(ShortTermMeanAcrossBirds(:,i)), std(ShortTermMeanAcrossBirds(:,i))/sqrt(size(ShortTermMeanAcrossBirds,1)), 'ks');
end
plot(repmat([1.1 1.9], size(ShortTermMeanAcrossBirds,1), 1)', ShortTermMeanAcrossBirds', 'ko-');
axis tight;
Temp = axis;
PlotAxis = [0.5 2.5 0 1.1*Temp(4)];
XAxis_LabelStrings = [{'Day 1'} {'1 day later'}];
LSINA_FixPlotLabels(gcf, 14, 'Fig_Position', [94 195 500 500], 'Fig_Title', '# of INs on consecutive days', 'YAxis_Label', 'Mean # of INs', 'Fig_XTickLabels', XAxis_LabelStrings, 'Fig_Axis', PlotAxis, 'Fig_Save', ['AllBirds.ShortTerm.NumINs.INBouts.png']);
p = signrank(ShortTermMeanAcrossBirds(:,1), ShortTermMeanAcrossBirds(:,2));
disp(['No of INs on Day 1 vs. No of INs on Day 2: Wilcoxon sign-rank test p is ', num2str(p)]);

% Plot summary across birds - long-term
figure;
hold on;
LongTermMeanAcrossBirds_Bar = bar(mean(LongTermMeanAcrossBirds));
set(LongTermMeanAcrossBirds_Bar, 'FaceColor', 'none', 'EdgeColor', 'k');
for i = 1:size(LongTermMeanAcrossBirds,2),
    errorbar(i, mean(LongTermMeanAcrossBirds(:,i)), std(LongTermMeanAcrossBirds(:,i))/sqrt(size(LongTermMeanAcrossBirds,1)), 'ks');
end
plot(repmat([1.1 1.9], size(LongTermMeanAcrossBirds,1), 1)', LongTermMeanAcrossBirds', 'ko-');
axis tight;
Temp = axis;
PlotAxis = [0.5 2.5 0 1.1*Temp(4)];
XAxis_LabelStrings = [{'Day 1'} {'6 months later'}];
LSINA_FixPlotLabels(gcf, 14, 'Fig_Position', [94 195 500 500], 'Fig_Title', '# of INs on days separated by 6 months', 'YAxis_Label', 'Mean # of INs', 'Fig_XTickLabels', XAxis_LabelStrings, 'Fig_Axis', PlotAxis, 'Fig_Save', ['AllBirds.LongTerm.NumINs.INBouts.png']);
p = signrank(ShortTermMeanAcrossBirds(:,1), ShortTermMeanAcrossBirds(:,2));
disp(['No of INs on Day 1 vs. No of INs 6 months later: Wilcoxon sign-rank test p is ', num2str(p)]);

% Now plot the acoustic properties of the first motif syllable depending on
% whether it was in an IN-only bout or other bouts

for i = 1:length(UniqueBirdIndices),
    FigLegend = [];
    figure;
    Indices = [];
    for j = 1:length(BirdParameters),
        if (~isempty(strfind(BirdParameters(j).BirdName, BirdParameters(UniqueBirdIndices(i)).BirdName)))
            Indices(end+1) = j;
            FigLegend{end+1} = BirdParameters(j).DataLabel;
        end
    end
    
    for j = 1:length(Indices),
        TotalSyllCount = 0;
        FirstMotifSyllIndex_INBouts = [];
        FirstMotifSyllIndex_NonINBouts = [];
        for k = 1:length(BirdParameters(Indices(j)).INResults.IndividualBoutLabels),
            FirstMotifSyll = find((BirdParameters(Indices(j)).INResults.AllCommonMotifPositions(1,:) > 0) & (BirdParameters(Indices(j)).INResults.AllCommonMotifPositions(2,:) == k), 1, 'first');
            FirstMotifSyll = FirstMotifSyll - TotalSyllCount;
            if ((FirstMotifSyll - 1) == sum(BirdParameters(Indices(j)).INResults.INs{k}(1:FirstMotifSyll-1)))
                FirstMotifSyllIndex_INBouts(end+1) = FirstMotifSyll + TotalSyllCount;
            else
                FirstMotifSyllIndex_NonINBouts(end+1) = FirstMotifSyll + TotalSyllCount;
            end
            TotalSyllCount = TotalSyllCount + length(BirdParameters(Indices(j)).INResults.IndividualBoutLabels{k});
        end
    end
end

% Now plot the number of INs for all birds as a summary with stats for two
% different types of bouts - ones beginning with INs and ones beginning
% with calls

% Now I want to find out the relationship between time to next IN and the
% change in acoustic properties. There are two possibilities
% 1) the time between INs matters - as time gets longer, then the next IN
% is more similar to the current IN
% 2) the time between INs does not matter - the next IN will always be at a
% fixed distance away from the current IN irrespective of the time to the
% next IN

% The way I'm going to do this is as follows:
% 1) First take all INs (irrespective of whether they're part of IN
% sequences before motifs or not
% 2) Then convert the feature space into principal components - I'm taking
% the 7 features listed below and will take the first 4 PCs
% 3) Then find repeated INs and measure time between them and the euclidean
% distance between them in the 4 PC space

% First convert everything into Principal Components
FeaturesToUse = [{'Duration'} {'LogAmplitude'} {'Entropy'} {'MeanFrequency'} {'PitchGoodness'} {'FrequencyModulation'}];

PCsToUse = [1 2 3 4];
RepeatNum = 10000;

for i = 1:length(BirdParameters),
    % Find all introductory notes
    AllINs = zeros(size(BirdParameters(i).AllLabels));
    for j = 1:length(BirdParameters(i).INLabels),
        AllINs(find(BirdParameters(i).AllLabels == BirdParameters(i).INLabels(j))) = 1;
    end
    AllINIndices = find(AllINs == 1);
    BirdParameters(i).AllINIndices = AllINIndices;
    
    AllSAPFeats = [];
    for j = 1:length(FeaturesToUse),
        AllSAPFeats(:,j) = eval(['BirdParameters(', num2str(i), ').SAPFeats.', FeaturesToUse{j}, '(:)']);
    end
    
    AllINSAPFeats = AllSAPFeats(AllINIndices,:);
    RawAllINSAPFeats = AllINSAPFeats;
    % Normalize by mean and standard deviation
    AllINSAPFeats = zscore(AllINSAPFeats);
    
    % Now get principal components
    [coeff, score, latent, tsquared, explained] = pca(AllINSAPFeats);
    
    disp(['Percent of total variance explained by first 4 PCs = ', num2str(sum(explained(PCsToUse))), '%']);
    BirdParameters(i).NeighbouringIN_RawFeatRatio = [];
    BirdParameters(i).NeighbouringIN_PCFeatRatio = [];
    BirdParameters(i).NeighbouringINs_TimesDistances = [];
    BirdParameters(i).NeighbouringIN_Indices = [];
    
    TotalSyllCount = 0;
    TotalINCount = 0;
    LastINs = find(BirdParameters(i).INResults.AllCommonMotifPreINs(end,AllINIndices) == -1);
    BirdParameters(i).LastIN_RawSAPFeats = RawAllINSAPFeats(LastINs,:);
    BirdParameters(i).MeanLastIN_RawSAPFeats = mean(RawAllINSAPFeats(LastINs,:));
    BirdParameters(i).STDLastIN_RawSAPFeats = std(RawAllINSAPFeats(LastINs,:));
    
    FirstINs = find(BirdParameters(i).INResults.AllCommonMotifPreINs(end-1,AllINIndices) == 1);
    BirdParameters(i).FirstIN_RawSAPFeats = RawAllINSAPFeats(FirstINs,:);
    BirdParameters(i).MeanFirstIN_RawSAPFeats = mean(RawAllINSAPFeats(FirstINs,:));
    BirdParameters(i).STDFirstIN_RawSAPFeats = std(RawAllINSAPFeats(FirstINs,:));
    
    BirdParameters(i).INPC_Scores = score;
    BirdParameters(i).MeanLastIN_PCScore = mean(score(LastINs, PCsToUse));
    BirdParameters(i).MeanFirstIN_PCScore = mean(score(FirstINs, PCsToUse));
    
    BirdParameters(i).FeaturesToUseForIN_PCs = FeaturesToUse;
    
    for j = 1:length(BirdParameters(i).INResults.INs),
        IN_Indices = find(BirdParameters(i).INResults.INs{j} == 1);
        for k = 1:length(IN_Indices)-1,
            if (IN_Indices(k+1) == (IN_Indices(k) + 1))
                BirdParameters(i).NeighbouringIN_Indices = [BirdParameters(i).NeighbouringIN_Indices; (TotalSyllCount + IN_Indices(k+1)); (TotalSyllCount + IN_Indices(k))];
                TimeBetweenINs = BirdParameters(i).AllOnsets(TotalSyllCount + IN_Indices(k+1)) - BirdParameters(i).AllOffsets(TotalSyllCount + IN_Indices(k));
                DistanceBetweenINs = pdist2(score(TotalINCount + k + 1, PCsToUse), score(TotalINCount + k, PCsToUse));
                DistanceToLastIN = [pdist2(score(TotalINCount + k, PCsToUse), BirdParameters(i).MeanLastIN_PCScore) pdist2(score(TotalINCount + k + 1, PCsToUse), BirdParameters(i).MeanLastIN_PCScore)];
                BirdParameters(i).NeighbouringINs_TimesDistances(end+1,:) = [TimeBetweenINs DistanceBetweenINs DistanceToLastIN BirdParameters(i).INResults.AllCommonMotifPreINs(end-1:end, TotalSyllCount + IN_Indices(k))'];
                BirdParameters(i).NeighbouringIN_RawFeatRatio(end+1,:) = RawAllINSAPFeats(TotalINCount + k + 1,:)./RawAllINSAPFeats(TotalINCount + k,:);
                BirdParameters(i).NeighbouringIN_PCFeatRatio(end+1,:) = score(TotalINCount + k + 1,PCsToUse)./score(TotalINCount + k,PCsToUse);
            end
        end
        TotalSyllCount = TotalSyllCount + length(BirdParameters(i).INResults.IndividualBoutLabels{j});
        TotalINCount = TotalINCount + length(IN_Indices);
    end
    BirdParameters(i).NeighbouringIN_Indices = unique(BirdParameters(i).NeighbouringIN_Indices);
end    

% Now to plot individual trials in PC space.
Colours = 'rbcmk';
Symbols = 'so^';

for i = 1:length(BirdParameters),
    figure;
    hold on;
    TempINPositions = BirdParameters(i).INResults.AllCommonMotifPreINs(end,:);
    PositionIndex = 0;
    for j = min(TempINPositions):1,
        Indices = find(TempINPositions == j);
        for k = 1:length(Indices),
            IN_Index = find(BirdParameters(i).AllINIndices == Indices(k));
            plot3(BirdParameters(i).INPC_Scores(IN_Index:(IN_Index + abs(j) - 1), 1), BirdParameters(i).INPC_Scores(IN_Index:(IN_Index + abs(j) - 1), 2), BirdParameters(i).INPC_Scores(IN_Index:(IN_Index + abs(j) - 1), 3), [Colours(mod(PositionIndex, length(Colours)) + 1), Symbols(mod(floor(PositionIndex/length(Colours)), length(Symbols)) + 1), '-']);;
            TempINPositions(Indices(k):(Indices(k) + abs(j) - 1)) = 0;
        end
        PositionIndex = PositionIndex + 1;
    end
end    
% Now I've got the change in ratio of feats. Now to simulate a simple
% model, where each IN starts at some state and then increases by the
% average ratio. If it is within the threshold, then the motif starts, else
% another IN is produced. Here I don't have any time in this model
for i = 1:length(BirdParameters),
    % Threshold is set as mean distance of all LastINs from the mean Last
    % IN acoustic structure + 2*std
    DistanceLastIN_FromMeanLastIN = pdist2(BirdParameters(i).LastIN_RawSAPFeats, mean(BirdParameters(i).LastIN_RawSAPFeats), 'mahalanobis', cov(BirdParameters(i).LastIN_RawSAPFeats));
    Threshold = mean(DistanceLastIN_FromMeanLastIN) + 2*std(DistanceLastIN_FromMeanLastIN);
    
    for j = 1:BirdParameters(i).NumBouts,
        % Get initial IN structure based on mean and std of first IN
        % acoustic structure
        Sim_IN_AcousticStructure = normrnd(BirdParameters(i).MeanFirstIN_RawSAPFeats, BirdParameters(i).STDFirstIN_RawSAPFeats);
        NumINs = 1;
        Flag = 0; % a flag which is 1 when IN structure is within 2 standard deviations of last IN structure and 0 otherwise
        % Now in a loop check if IN structure is within 2 standard
        % deviations of Last IN structure
        while (Flag == 0)
            MahalDistanceToLastIN = pdist2(Sim_IN_AcousticStructure, BirdParameters(i).MeanLastIN_RawSAPFeats, 'mahalanobis', cov(BirdParameters(i).LastIN_RawSAPFeats));
            if (MahalDistanceToLastIN <= Threshold)
                Flag = 1;
            else
                Increment = normrnd(mean(BirdParameters(i).NeighbouringIN_RawFeatRatio), std(BirdParameters(i).NeighbouringIN_RawFeatRatio));
                Sim_IN_AcousticStructure = Sim_IN_AcousticStructure.*Increment;
                NumINs = NumINs + 1;
            end
        end
        BirdParameters(i).AcousticModel_SimNumINs(j) = NumINs;
    end
end
% Now get parameters for a model. Basically get the following parameters
% 1) Initial state in terms of acoustic parameters of IN and interval to
% next IN
% 2) Relative change in acoustic state and interval as the bird sings
% another IN
% 3) Threshold in terms of acoustic parameters and interval to motif to
% determine when IN sequence is going to end

% Use only 4 acoustic features - duration, mean frequency, log amplitude,
% entropy
FeaturesToUse = [{'Duration'} {'LogAmplitude'} {'Entropy'} {'MeanFrequency'}];

for i = 1:length(BirdParameters),
    BirdParameters(i).InitialInterval = [];
    BirdParameters(i).InitialAcousticState = [];
    
    BirdParameters(i).IntervalRatio = [];
    BirdParameters(i).AcousticStateRatio = [];
    
    BirdParameters(i).LastInterval = [];
    BirdParameters(i).LastAcousticState = [];
    
    for j = 1:length(BirdParameters(i).IN_NonIN_BoutIndices),
        if (BirdParameters(i).IN_NonIN_BoutIndices(j) == 1)
            InitialIN_Index = find((BirdParameters(i).INResults.AllBoutStartINs(2,:) == j) & (BirdParameters(i).INResults.AllBoutStartINs(3,:) == 1));
            BirdParameters(i).InitialInterval(end+1) = BirdParameters(i).AllOnsets(InitialIN_Index +1) - BirdParameters(i).AllOffsets(InitalIN_Index);
            TempSAPFeats = [];
            for k = 1:length(FeaturesToUse),
                TempSAPFeats(1,k) = eval(['BirdParameters(', num2str(i), ').SAPFeats.', FeaturesToUse{k}, '(', num2str(InitialIN_Index), ')']);
            end
            BirdParameters(i).InitialAcousticState(end+1,:) = TempSAPFeats;
            
            LastIN_Index = find((BirdParameters(i).INResults.AllBoutStartINs(2,:) == j) & (BirdParameters(i).INResults.AllBoutStartINs(4,:) == -1));
            BirdParameters(i).LastInterval(end+1) = BirdParameters(i).AllOnsets(LastIN_Index +1) - BirdParameters(i).AllOffsets(LastIN_Index);
            TempSAPFeats = [];
            for k = 1:length(FeaturesToUse),
                TempSAPFeats(1,k) = eval(['BirdParameters(', num2str(i), ').SAPFeats.', FeaturesToUse{k}, '(', num2str(LastIN_Index), ')']);
            end
            BirdParameters(i).LastAcousticState(end+1,:) = TempSAPFeats;
            
            for NumINs = 1:BirdParameters(i).INResults.NumBoutStartINs,
            end
        end
    end
end
    
for i = 1:length(UniqueBirdIndices),
end
disp('Finished Analysis');