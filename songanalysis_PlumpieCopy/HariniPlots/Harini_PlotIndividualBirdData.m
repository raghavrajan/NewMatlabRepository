function [BirdStats] = Harini_PlotIndividualBirdData(UniqueBirds, BirdParameters, Date, DateNumber, Condition)

BirdStats.BirdParameters = BirdParameters;
BirdStats.BirdName = UniqueBirds;
BirdStats.Date = Date;
BirdStats.Condition = Condition;
BirdStats.DateNumber = DateNumber;

% First plot the number of valid bouts for each day along with the number
% and proportion of D, DUN and UN bouts for that day.

figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [5 494 1850 500]);
p = panel();
p.pack({1});

BoutTypes = {'D' 'DUN' 'UN'};
ResponseTypes = {'R' 'NR'};
XLabelString = [];
for i = 1:length(Date),
    NumBouts(i) = length(BirdParameters(i).ValidSongBouts);
    BirdStats.LatencyToFirstSong(i) = BirdParameters(i).LatencyToFirstSong;
    BirdStats.LatencyToFirstMotif(i) = BirdParameters(i).LatencyToFirstMotif;
    for j = 1:length(BoutTypes),
        PropBoutTypes(i,j) = length(strmatch(BoutTypes{j}, BirdParameters(i).BoutVideoScoringTag, 'exact'));
    end
    for j = 1:length(ResponseTypes),
        PropResponseTypes(i,j) = length(strmatch(ResponseTypes{j}, BirdParameters(i).BoutFemaleResponseTag, 'exact'));
    end
    XLabelString{i} = [Date{i}, ' : ', Condition{i}];
end

p(1).select();
hold on;
bar((1:1:length(Date)) - 0.25, NumBouts, 0.2)
bar((1:1:length(Date)), PropBoutTypes, 0.2, 'stacked');
bar((1:1:length(Date)) + 0.25, PropResponseTypes, 0.2, 'stacked');
axis tight;
Temp = axis;
Temp = [0.5 (length(Date)+0.5) 0 1.2*Temp(4)];
axis(Temp);
set(gca, 'XTick', 1:1:length(Date), 'XTickLabel', XLabelString, 'XTickLabelRotation', 30);
xlabel('Recording Day');
ylabel('Number of bouts');
legend('Total valid bouts', 'D bouts', 'DUN bouts', 'UN bouts', 'R bouts', 'NR bouts');
p.fontsize = 16;
p.marginbottom = 35;
p.marginleft = 20;

FinalFigureDir = '/home/raghav/StudentRelated/Harini/Manuscript/JEBManuscript/Figures/';

print(fullfile(FinalFigureDir, [BirdStats.BirdName, '.AllDaysBoutNumbers.png']), '-dpng', '-r300');


% Now to get all stats
% First # of days
BirdStats.NumRecordingDays = length(unique(Date));
BirdStats.TotalNumDaysFromFirstToLast = max(DateNumber) - min(DateNumber) + 1;
BirdStats.UniqueConditions = unique(Condition);
BirdStats.NumUniqueConditions = length(unique(Condition));
for i = 1:length(BirdStats.UniqueConditions),
    BirdStats.NumSessionsPerCondition(i) = length(strmatch(BirdStats.UniqueConditions{i}, Condition, 'exact'));
end
BirdStats.TotalNumSessions = sum(BirdStats.NumSessionsPerCondition);
BirdStats.TotalNumSongBouts = sum(NumBouts);
% Now to find the % of bouts that are continuous with next - these won't be
% used
for i = 1:length(Date),
    BirdStats.NumAllSongBouts(i) = length(BirdParameters(i).BoutContinuousWithNext);
    BirdStats.NumContinuousSongBouts(i) = length(find(BirdParameters(i).BoutContinuousWithNext == 1));
end
% Fraction of bouts thrown away
BirdStats.PercentContinousSongBouts = 100 * sum(BirdStats.NumContinuousSongBouts)/sum(BirdStats.NumAllSongBouts);

% Now to find out fraction of bouts with video scoring after excluding the
% continuous song bouts
for i = 1:length(Date),
    SongBoutsWithNoVideoScoring = strmatch('NA', BirdParameters(i).BoutVideoScoringTag, 'exact');
    SongBoutsWithVideoScoring = setdiff((1:1:length(BirdParameters(i).BoutVideoScoringTag)), SongBoutsWithNoVideoScoring);
    NonContinuousSongBouts = find(BirdParameters(i).BoutContinuousWithNext == 0);
    BirdStats.NumNonContinuousSongBoutsWithVideoScoring(i) = length(intersect(NonContinuousSongBouts, SongBoutsWithVideoScoring));
    if (isempty(strmatch('UN', Condition{i}, 'exact')))
        % Now to find all of the 'D', 'DUN' and 'UN' non-continuous song bouts
        % each session
        for j = 1:length(BoutTypes),
            BirdStats.NumBoutTypes(i,j) = length(intersect(strmatch(BoutTypes{j}, BirdParameters(i).BoutVideoScoringTag, 'exact'), NonContinuousSongBouts));
        end
    else
        BirdStats.NumBoutTypes(i,3) = length(NonContinuousSongBouts);
    end
end

for i = 1:length(BirdStats.UniqueConditions),
    Indices = strmatch(BirdStats.UniqueConditions{i}, Condition, 'exact');
    BirdStats.NumBoutTypesPerCondition(i,:) = sum(BirdStats.NumBoutTypes(Indices,:));
    BirdStats.LatencyToFirstSongPerCondition{i} = BirdStats.LatencyToFirstSong(Indices);
    BirdStats.LatencyToFirstMotifPerCondition{i} = BirdStats.LatencyToFirstMotif(Indices);
end



disp('Finished plotting');
