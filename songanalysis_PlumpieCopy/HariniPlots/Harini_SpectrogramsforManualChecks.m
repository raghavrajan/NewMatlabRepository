function [] = Harini_SpectrogramsforManualChecks(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)

% =============== Plots for manual checks =================================
% THis script is to generate some plots for manually verifying that the
% program is indeed calcuating things properly. So, what I want to do is
% plot out spectrograms for each of the sessions with the labels and times
% of the individual syllables - different colours for motif syllables and
% different colours for INs and another colour for all other syllables
% Then I need to mark in red the bouts that have been used and in blue the
% bouts that haven't been used.
% Then mark the IN numbers and complete and partial motif numbers.
% =========================================================================
FigureDir = '/data3/raghav/Harini_Plots';
FontSizeVal = 10;
NumRowsPerFigure = 6;

[DirBout_Stats, UndirBout_Stats, DirUndirBout_Stats] = Harini_ProcessINNos(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones);

for i = 1:length(IndividualBirds),
    BoutIndexColumnIndex = find(cellfun(@length,strfind(IndividualBirds(i).BoutStatisticsColumnNames, 'BoutIndex')));
    NumINsColumnIndex = find(cellfun(@length,strfind(IndividualBirds(i).BoutStatisticsColumnNames, 'NumINs')));
    NumMotifsColumnIndex = find(cellfun(@length,strfind(IndividualBirds(i).BoutStatisticsColumnNames, 'NumMotifs')));
    FMDColumnIndex = find(cellfun(@length,strfind(IndividualBirds(i).BoutStatisticsColumnNames, 'FirstMotifDur')));
    
    Index = 1;
    
    INLabels = IndividualBirds(i).SortedBirdParameters(1).INLabels;
    for j = 1:length(INLabels),
        INLabelCase{j} = INLabels(j);
    end
    
    MotifLabels = IndividualBirds(i).SortedBirdParameters(1).MotifLabels;
    for j = 1:length(MotifLabels),
        MotifLabelCase{j} = MotifLabels(j);
    end
    
    CommonMotifs = IndividualBirds(i).SortedBirdParameters(1).CommonMotifs;
    
    ValidSongBouts = find((IndividualBirds(i).Bouts(:,7) == 1) & (IndividualBirds(i).Bouts(:,8) > 0) & (IndividualBirds(i).Bouts(:,9) > 1));
    for j = 1:length(IndividualBirds(i).SongFileNames),
        if (mod(Index, NumRowsPerFigure) == 1)
            if (Index > 1)
                p.fontsize = FontSizeVal;
    %             p.margintop = 10;
    %             p.marginleft = 20;
    %             p.marginright = 20;

                set(gcf, 'PaperPositionMode', 'auto');
                print(fullfile(FigureDir, [BirdNames{i}, '.Spectrograms.#', num2str(Index), '.png']), '-dpng', '-r300');
                close(gcf);
            end
            
            figure;
            set(gcf, 'Position', [1 33 1920 964]);
            set(gcf, 'Color', 'w');
            p = panel();
            clear FigureRows;
            for Rows = 1:NumRowsPerFigure,
                FigureRows{Rows} = 1/NumRowsPerFigure;
            end
            p.pack(FigureRows);
            p.de.margin = 4;
        end
        [FileDir, FileName, FileExt] = fileparts(IndividualBirds(i).SongFileNames{j});
        [Data, Fs] = GetData(FileDir, [FileName, FileExt], IndividualBirds(1).SortedBirdParameters(1).FileType, 0);
        Bouts = find(IndividualBirds(i).Bouts(:,3) == j);
        if (~isempty(Bouts))
            BoutOnset = IndividualBirds(i).Bouts(Bouts,5);
            BoutOffset = IndividualBirds(i).Bouts(Bouts,6);
        
            BoutSyllLabels = char(IndividualBirds(i).AllSyllableData(IndividualBirds(i).Bouts(Bouts(1),1):IndividualBirds(i).Bouts(Bouts(end),2),1));
            BoutSyllOnsets = IndividualBirds(i).AllSyllableData(IndividualBirds(i).Bouts(Bouts(1),1):IndividualBirds(i).Bouts(Bouts(end),2),4);
            BoutSyllOffsets = IndividualBirds(i).AllSyllableData(IndividualBirds(i).Bouts(Bouts(1),1):IndividualBirds(i).Bouts(Bouts(end),2),5);
            
            BoutDistance = mean(IndividualBirds(i).AllConditionIndices(IndividualBirds(i).Bouts(Bouts(1),1):IndividualBirds(i).Bouts(Bouts(end),2)));
            BoutFirstDistance = IndividualBirds(i).AllConditionIndices(IndividualBirds(i).Bouts(Bouts(1),1));
            switch BoutFirstDistance
                case 1
                    BoutFirstDistanceString = 'L0';
                    
                case 2
                    BoutFirstDistanceString = 'L1';
                    
                case 3
                    BoutFirstDistanceString = 'L2';
                    
                case 4
                    BoutFirstDistanceString = 'L3';
                    
                case 5
                    BoutFirstDistanceString = 'L4';
                    
                case 6
                    BoutFirstDistanceString = 'UN';
            end
            
            % Now get the # of INs as calculated by bout statistics and by
            % the other script
            
        end
        for k = 1:3,
            switch (k)
                case 1
                    PlotData = Data(1:round(length(Data)/3));
                    Time = (1:1:length(PlotData))/Fs;
                case 2
                    PlotData = Data((round(length(Data)/3) + 1):(2*(round(length(Data)/3))));
                    Time = Time(end) + (1:1:length(PlotData))/Fs;
                case 3
                    PlotData = Data(((2*round(length(Data)/3)) + 1):end);
                    Time = Time(end) + (1:1:length(PlotData))/Fs;
            end
        
            if (mod(Index, NumRowsPerFigure) == 0)
                p(NumRowsPerFigure).select();
            else
                p(Index - floor(Index/NumRowsPerFigure)*NumRowsPerFigure).select();
            end
            PlotSpectrogramInAxis_SongVar(PlotData, Time, Fs, gca);
            set(gca, 'Visible', 'off');
            TempFileName = FileName;
            TempFileName(find(FileName == '_')) = ' ';
            text(Time(1), 500, [TempFileName, FileExt], 'FontSize', 10);
            
            for BoutNo = 1:length(Bouts),
                if (~isempty(intersect(ValidSongBouts, Bouts(BoutNo))))
                    PlotColour = 'r';
                    BoutStatisticsBoutNo = find(IndividualBirds(i).BoutStatistics(:, BoutIndexColumnIndex) == Bouts(BoutNo));
                    if (~isempty(BoutStatisticsBoutNo))
                        BoutType = IndividualBirds(i).BoutCategorisation{BoutStatisticsBoutNo};
                        BoutStatsNumINs = IndividualBirds(i).BoutStatistics(BoutStatisticsBoutNo, NumINsColumnIndex);
                        BoutStatsNumMotifs = IndividualBirds(i).BoutStatistics(BoutStatisticsBoutNo, NumMotifsColumnIndex);
                        BoutStatsFMD = IndividualBirds(i).BoutStatistics(BoutStatisticsBoutNo, FMDColumnIndex(1));
                        BoutStatsString = ['BS: #INs=', num2str(BoutStatsNumINs), ', #Motifs=', num2str(BoutStatsNumMotifs), ', FMD=', num2str(BoutStatsFMD), 'ms'];
                    else
                        BoutStatsString = ['BS: Bout not found'];
                    end
                    FoundBoutFlag = 0;
                    for DirBoutNo = 1:length(DirBout_Stats{i}),
                        if (~isempty(DirBout_Stats{i}{DirBoutNo}))
                            if (~isempty(find([DirBout_Stats{i}{DirBoutNo}.BoutIndex] == Bouts(BoutNo))))
                                FoundBoutFlag = 1;
                                FoundBoutIndex = find([DirBout_Stats{i}{DirBoutNo}.BoutIndex] == Bouts(BoutNo));
                                OtherStatsNumINs = DirBout_Stats{i}{DirBoutNo}.TotalINNumber_500ms(FoundBoutIndex);
                                OtherStatsNumINs = DirBout_Stats{i}{DirBoutNo}.TotalINNumber_500ms(FoundBoutIndex);
                                OtherStatsNumFullMotifs = DirBout_Stats{i}{DirBoutNo}.CompleteMotifNumber(FoundBoutIndex);
                                OtherStatsNumPartialMotifs = DirBout_Stats{i}{DirBoutNo}.PartialMotifNumber(FoundBoutIndex);
                                OtherStatsFMD = DirBout_Stats{i}{DirBoutNo}.FirstMotifDuration(FoundBoutIndex);
                                OtherStatsString = ['OS: D / #INs=', num2str(OtherStatsNumINs), ', #FullMotifs=', num2str(OtherStatsNumFullMotifs), ', #PartialMotifs=', num2str(OtherStatsNumPartialMotifs), ', FMD=', num2str(OtherStatsFMD), 'ms'];
                                break;
                            end
                        end
                    end
                    
                    if (FoundBoutFlag == 0)
                        for DirUndirBoutNo = 1:length(DirUndirBout_Stats{i}),
                            if (~isempty(DirUndirBout_Stats{i}{DirUndirBoutNo}))
                                if (~isempty(find([DirUndirBout_Stats{i}{DirUndirBoutNo}.BoutIndex] == Bouts(BoutNo))))
                                    FoundBoutFlag = 1;
                                    FoundBoutIndex = find([DirUndirBout_Stats{i}{DirUndirBoutNo}.BoutIndex] == Bouts(BoutNo));
                                    OtherStatsNumINs = DirUndirBout_Stats{i}{DirUndirBoutNo}.TotalINNumber_500ms(FoundBoutIndex);
                                    OtherStatsNumINs = DirUndirBout_Stats{i}{DirUndirBoutNo}.TotalINNumber_500ms(FoundBoutIndex);
                                    OtherStatsNumFullMotifs = DirUndirBout_Stats{i}{DirUndirBoutNo}.CompleteMotifNumber(FoundBoutIndex);
                                    OtherStatsNumPartialMotifs = DirUndirBout_Stats{i}{DirUndirBoutNo}.PartialMotifNumber(FoundBoutIndex);
                                    OtherStatsFMD = DirUndirBout_Stats{i}{DirUndirBoutNo}.FirstMotifDuration(FoundBoutIndex);
                                    OtherStatsString = ['OS: DUN / #INs=', num2str(OtherStatsNumINs), ', #FullMotifs=', num2str(OtherStatsNumFullMotifs), ', #PartialMotifs=', num2str(OtherStatsNumPartialMotifs), ', FMD=', num2str(OtherStatsFMD), 'ms'];
                                    break;
                                end
                            end
                        end
                    end
                    
                    if (FoundBoutFlag == 0)
                        for UndirBoutNo = 1:length(UndirBout_Stats{i}),
                            if (~isempty(UndirBout_Stats{i}{UndirBoutNo}))
                                if (~isempty(find([UndirBout_Stats{i}{UndirBoutNo}.BoutIndex] == Bouts(BoutNo))))
                                    FoundBoutFlag = 1;
                                    FoundBoutIndex = find([UndirBout_Stats{i}{UndirBoutNo}.BoutIndex] == Bouts(BoutNo));
                                    OtherStatsNumINs = UndirBout_Stats{i}{UndirBoutNo}.TotalINNumber_500ms(FoundBoutIndex);
                                    OtherStatsNumINs = UndirBout_Stats{i}{UndirBoutNo}.TotalINNumber_500ms(FoundBoutIndex);
                                    OtherStatsNumFullMotifs = UndirBout_Stats{i}{UndirBoutNo}.CompleteMotifNumber(FoundBoutIndex);
                                    OtherStatsNumPartialMotifs = UndirBout_Stats{i}{UndirBoutNo}.PartialMotifNumber(FoundBoutIndex);
                                    OtherStatsFMD = UndirBout_Stats{i}{UndirBoutNo}.FirstMotifDuration(FoundBoutIndex);
                                    OtherStatsString = ['OS: UN / #INs=', num2str(OtherStatsNumINs), ', #FullMotifs=', num2str(OtherStatsNumFullMotifs), ', #PartialMotifs=', num2str(OtherStatsNumPartialMotifs), ', FMD=', num2str(OtherStatsFMD), 'ms'];
                                    break;
                                end
                            end
                        end
                    end
                    if (FoundBoutFlag == 0)
                        OtherStatsString = ['OS: No bout found'];
                    end
                else
                    PlotColour = 'b';
                    BoutType = ' ';
                    BoutStatsString = ' ';
                    OtherStatsString = ' ';
                end
                if ((BoutOnset(BoutNo)/1000 >= Time(1)) && (BoutOnset(BoutNo)/1000 <= Time(end)))
                    if (BoutOffset(BoutNo)/1000 <= Time(end))
                        plot([BoutOnset(BoutNo) BoutOffset(BoutNo)]/1000, [7500 7500], PlotColour, 'LineWidth', 1);
                        text(BoutOnset(BoutNo)/1000, 2000+(floor(BoutNo/3)*1000), [num2str(Bouts(BoutNo)), ':', BoutType, ':', num2str(BoutDistance), ':', BoutFirstDistanceString, ':', BoutStatsString, ':', OtherStatsString], 'Color', 'r', 'FontSize', 8, 'FontWeight', 'bold');
                    else
                        plot([BoutOnset(BoutNo) Time(end)*1000]/1000, [7500 7500], PlotColour, 'LineWidth', 1);
                        if (BoutOnset(BoutNo)/1000 <= (Time(1) + ((Time(end) - Time(1))/2)))
                            text(BoutOnset(BoutNo)/1000, 2000+(floor(BoutNo/3)*1000), [num2str(Bouts(BoutNo)), ':', BoutType, ':', num2str(BoutDistance), ':', BoutFirstDistanceString, ':', BoutStatsString, ':', OtherStatsString], 'Color', 'r', 'FontSize', 8, 'FontWeight', 'bold');
                        else
                            text(Time(1) + ((Time(end)-Time(1))/2), 2000+(floor(BoutNo/3)*1000), [num2str(Bouts(BoutNo)), ':', BoutType, ':', num2str(BoutDistance), ':', BoutFirstDistanceString, ':', BoutStatsString, ':', OtherStatsString], 'Color', 'r', 'FontSize', 8, 'FontWeight', 'bold');
                        end
                    end
                else
                    if ((BoutOffset(BoutNo)/1000 >= Time(1)) && (BoutOffset(BoutNo)/1000 <= Time(end)))
                        plot([Time(1)*1000 BoutOffset(BoutNo)]/1000, [7500 7500], PlotColour, 'LineWidth', 1);
                    end
                end
            end
            if (~isempty(BoutSyllLabels))
                for SyllNo = 1:length(BoutSyllLabels),
                    switch (BoutSyllLabels(SyllNo))
                        case INLabelCase
                            TextColor = 'r';
                            BoutOnsetText = num2str(round(BoutSyllOnsets(SyllNo)));
                        case MotifLabelCase
                            TextColor = 'b';
                            BoutOnsetText = '';
                        otherwise
                            TextColor = 'c';
                            BoutOnsetText = '';
                    end
                    
                    if ((BoutSyllOnsets(SyllNo)/1000 >= Time(1)) && (BoutSyllOnsets(SyllNo)/1000 <= Time(end)))
                        text(BoutSyllOnsets(SyllNo)/1000, 6500, [BoutSyllLabels(SyllNo)], 'FontSize', 15, 'FontWeight', 'bold', 'Color', TextColor);
                        text(BoutSyllOnsets(SyllNo)/1000, 5500, BoutOnsetText, 'FontSize', 6, 'FontWeight', 'bold', 'Color', TextColor);
                    end
                end
            end
            Index = Index + 1;
        end
    end
end
        
disp('Done');