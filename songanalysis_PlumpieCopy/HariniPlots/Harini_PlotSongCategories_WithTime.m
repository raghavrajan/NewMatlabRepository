function [DirBout_Stats, UndirBout_Stats, DirUndirBout_Stats] = Harini_PlotSongCategories_WithTime(IndividualBirds, BirdNames, BirdOption, FeatureToPlot_ColNo, PlotType, DifferentMicrophones)


FinalFigureDir = '/home/raghav/StudentRelated/Harini/Manuscript/results';

TimeBinSize = 5; % in minutes

for i = 1:length(IndividualBirds),
    ConditionColumnIndex = find(cellfun(@length,strfind(IndividualBirds(i).BoutStatisticsColumnNames, 'Condition')));
    BoutIndexColumnIndex = find(cellfun(@length,strfind(IndividualBirds(i).BoutStatisticsColumnNames, 'BoutIndex')));
    figure;
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [463 41 750 950]);
    p = panel();
    p.pack(6,1);
    p.de.margin = 5;
    
    for j = min(IndividualBirds(i).BoutStatistics(:,ConditionColumnIndex)):max(IndividualBirds(i).BoutStatistics(:,ConditionColumnIndex)),
        AllBouts{i}{j} = [];
        Bouts{i}{j} = IndividualBirds(i).BoutStatistics(find(IndividualBirds(i).BoutStatistics(:,ConditionColumnIndex) == j), BoutIndexColumnIndex);
        
        DirBouts{i}{j} = intersect(IndividualBirds(i).BoutStatistics(strmatch('D', IndividualBirds(i).BoutCategorisation, 'exact'), BoutIndexColumnIndex), Bouts{i}{j});
        if (~isempty(DirBouts{i}{j}))
            AllBouts{i}{j} = [AllBouts{i}{j}; [DirBouts{i}{j}(:) ones(size(DirBouts{i}{j}(:)))*1]];
        end
        
        DirUndirBouts{i}{j} = intersect(IndividualBirds(i).BoutStatistics(strmatch('DUN', IndividualBirds(i).BoutCategorisation, 'exact'), BoutIndexColumnIndex), Bouts{i}{j});
        if (~isempty(DirUndirBouts{i}{j}))
            AllBouts{i}{j} = [AllBouts{i}{j}; [DirUndirBouts{i}{j}(:) ones(size(DirUndirBouts{i}{j}(:)))*2]];
        end
        
        if (j ~= 6)
            UndirBouts{i}{j} = intersect(IndividualBirds(i).BoutStatistics(strmatch('UN', IndividualBirds(i).BoutCategorisation, 'exact'), BoutIndexColumnIndex), Bouts{i}{j});
        else
            UndirBouts{i}{j} = Bouts{i}{j}; % since these bouts are anyway only in the absence of the female.
        end
        if (~isempty(UndirBouts{i}{j}))
            AllBouts{i}{j} = [AllBouts{i}{j}; [UndirBouts{i}{j}(:) ones(size(UndirBouts{i}{j}(:)))*3]];
        end
        
        p(j,1).select();
        if (~isempty(AllBouts{i}{j}))
            AllBoutTimes{i}{j} = IndividualBirds(i).BoutTimesRelativeToStart(AllBouts{i}{j}(:,1)) * 60; % in minutes
            Index = 1;
            for Time = 0:TimeBinSize:120,
                Songs = find((AllBoutTimes{i}{j} > Time) & (AllBoutTimes{i}{j} <= (Time + TimeBinSize)));
                if (~isempty(Songs))
                    for SongType = 1:3,
                        PropSongs{i}{j}(Index,SongType) = length(find(AllBouts{i}{j}(Songs,2) == SongType))/length(Songs);
                    end
                else
                    PropSongs{i}{j}(Index,1:3) = 0;
                end
                Index = Index + 1;
            end
            bar((0:TimeBinSize:120) + TimeBinSize/2, PropSongs{i}{j}, 'stacked');
            legend('D', 'DUN', 'UN');
        end
        
        axis([0 140 0 1]);
        if (j == 1)
            title([BirdNames{i}, ': Proportion of songs scored as D, DUN and UN at different times after start']);
        end
        if (j == 6)
            xlabel('Time bins (minutes)');
        else
            set(gca, 'XTickLabel', []);
        end
        if (j == 3)
            ylabel('Fraction of songs');
        end
        set(gca, 'Box', 'off');
    end
    p.fontsize = 12;
    p.margintop = 10;
    set(gcf, 'ReSize', 'off');
    set(gcf, 'PaperPositionMode', 'auto');
    print(fullfile(FinalFigureDir, [BirdNames{i}, '.ProportionofSongTypesRelativeToStartTime.eps']), '-depsc2', '-r600');
    print(fullfile(FinalFigureDir, [BirdNames{i}, '.ProportionofSongTypesRelativeToStartTime.png']), '-dpng', '-r600');
end
disp('Finished');

