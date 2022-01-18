function [BirdStats] = Harini_GetIndividualBirdNumINs_Motifs_All(BirdStats, InterBoutInterval, BirdId)

BoutTypes = {'D' 'DUN' 'UN'};
ResponseTypes = {'R' 'NR'};
% Now I need to find all the valid bouts (non continuous with video
% scoring) and get all the parameters for those bouts - # of INs, # of
% motifs, etc.

BirdStats.BirdId = BirdId;
BirdStats.RecordingDayIndex = dummyvar((BirdStats.DateNumber - BirdStats.DateNumber(1))+1) * (1:1:(BirdStats.DateNumber(end) - BirdStats.DateNumber(1) + 1))';

for i = 1:length(BirdStats.Condition),
    switch (BirdStats.Condition{i})
        case 'L0'
            BirdStats.ConditionIndex(i) = 1;
        case 'L1'
            BirdStats.ConditionIndex(i) = 2;
        case 'L2'
            BirdStats.ConditionIndex(i) = 3;
        case 'L3'
            BirdStats.ConditionIndex(i) = 4;
        case 'L4'
            BirdStats.ConditionIndex(i) = 5;
        case 'UN'
            BirdStats.ConditionIndex(i) = 6;
    end
end

for i = 1:length(BirdStats.Date),
    NonContinuousSongBouts = find((BirdStats.BirdParameters(i).BoutContinuousWithNext == 0) & (BirdStats.BirdParameters(i).BoutPreTime >= InterBoutInterval/1000) & (BirdStats.BirdParameters(i).BoutPostTime >= InterBoutInterval/1000));
    if (isempty(strmatch('UN', BirdStats.Condition{i}, 'exact')))
        % Now to find all of the 'D', 'DUN' and 'UN' non-continuous song bouts
        % each session
        for j = 1:length(BoutTypes),
            switch (BoutTypes{j})
                case 'D'
                    BirdStats.D_songs{i} = (intersect(strmatch(BoutTypes{j}, BirdStats.BirdParameters(i).BoutVideoScoringTag, 'exact'), NonContinuousSongBouts));
                    
                case 'DUN'
                    BirdStats.DUN_songs{i} = (intersect(strmatch(BoutTypes{j}, BirdStats.BirdParameters(i).BoutVideoScoringTag, 'exact'), NonContinuousSongBouts));
                    
                case 'UN'
                    BirdStats.UN_songs{i} = (intersect(strmatch(BoutTypes{j}, BirdStats.BirdParameters(i).BoutVideoScoringTag, 'exact'), NonContinuousSongBouts));
            end
        end
    else
        BirdStats.UN_songs{i} = NonContinuousSongBouts;
    end
end

% Now with all the bouts chosen, should get IN numbers - will pool all of
% them for now
% Also, have a variable that has all of the data for Num INs in the form of
% a big matrix - first column will be number of INs, the next columns will
% be various group descriptors - distance first, then session #, then
% recording day #, then time of day, bird identity
BirdStats.AllValidNumINs = [];
BirdStats.AllValidNumMotifs = [];
BirdStats.AllValidFirstMotifDur = [];

for i = 1:length(BirdStats.UniqueConditions),
    Indices = strmatch(BirdStats.UniqueConditions{i}, BirdStats.Condition, 'exact');
    if (isempty(strmatch(BirdStats.UniqueConditions{i}, 'UN', 'exact')))
        BirdStats.ConditionSpecific_DSong_INs{i} = [];
        BirdStats.ConditionSpecific_DSong_NumMotifs{i} = [];
        BirdStats.ConditionSpecific_DSong_FirstMotifDur{i} = [];
        
        BirdStats.ConditionSpecific_DSong_MeanINs{i} = [];
        BirdStats.ConditionSpecific_DSong_MeanNumMotifs{i} = [];
        BirdStats.ConditionSpecific_DSong_MeanFirstMotifDur{i} = [];

        for j = 1:length(Indices),
            BirdStats.ConditionSpecific_DSong_INs{i} = [BirdStats.ConditionSpecific_DSong_INs{i} BirdStats.BirdParameters(Indices(j)).NumINs(BirdStats.D_songs{Indices(j)})];
            BirdStats.AllValidNumINs = [BirdStats.AllValidNumINs; [BirdStats.BirdParameters(Indices(j)).NumINs(BirdStats.D_songs{Indices(j)})' ones(length(BirdStats.D_songs{Indices(j)}),1)*BirdStats.ConditionIndex(Indices(j)) ones(length(BirdStats.D_songs{Indices(j)}),1)*Indices(j) ones(length(BirdStats.D_songs{Indices(j)}),1)*BirdStats.RecordingDayIndex(Indices(j)) ones(length(BirdStats.D_songs{Indices(j)}),1)*BirdStats.BirdId]];
            
            BirdStats.ConditionSpecific_DSong_NumMotifs{i} = [BirdStats.ConditionSpecific_DSong_NumMotifs{i} BirdStats.BirdParameters(Indices(j)).NumMotifs(BirdStats.D_songs{Indices(j)})];
            BirdStats.AllValidNumMotifs = [BirdStats.AllValidNumMotifs; [BirdStats.BirdParameters(Indices(j)).NumMotifs(BirdStats.D_songs{Indices(j)})' ones(length(BirdStats.D_songs{Indices(j)}),1)*BirdStats.ConditionIndex(Indices(j)) ones(length(BirdStats.D_songs{Indices(j)}),1)*Indices(j) ones(length(BirdStats.D_songs{Indices(j)}),1)*BirdStats.RecordingDayIndex(Indices(j)) ones(length(BirdStats.D_songs{Indices(j)}),1)*BirdStats.BirdId]];
            
            BirdStats.ConditionSpecific_DSong_FirstMotifDur{i} = [BirdStats.ConditionSpecific_DSong_FirstMotifDur{i} BirdStats.BirdParameters(Indices(j)).FirstMotifDur(BirdStats.D_songs{Indices(j)})];
            BirdStats.AllValidFirstMotifDur = [BirdStats.AllValidFirstMotifDur; [BirdStats.BirdParameters(Indices(j)).FirstMotifDur(BirdStats.D_songs{Indices(j)})' ones(length(BirdStats.D_songs{Indices(j)}),1)*BirdStats.ConditionIndex(Indices(j)) ones(length(BirdStats.D_songs{Indices(j)}),1)*Indices(j) ones(length(BirdStats.D_songs{Indices(j)}),1)*BirdStats.RecordingDayIndex(Indices(j)) ones(length(BirdStats.D_songs{Indices(j)}),1)*BirdStats.BirdId]];
            
            BirdStats.ConditionSpecific_DSong_MeanINs{i} = [BirdStats.ConditionSpecific_DSong_MeanINs{i} mean(BirdStats.BirdParameters(Indices(j)).NumINs(BirdStats.D_songs{Indices(j)}))];
            BirdStats.ConditionSpecific_DSong_MeanNumMotifs{i} = [BirdStats.ConditionSpecific_DSong_MeanNumMotifs{i} mean(BirdStats.BirdParameters(Indices(j)).NumMotifs(BirdStats.D_songs{Indices(j)}))];
            BirdStats.ConditionSpecific_DSong_MeanFirstMotifDur{i} = [BirdStats.ConditionSpecific_DSong_MeanFirstMotifDur{i} nanmean(BirdStats.BirdParameters(Indices(j)).FirstMotifDur(BirdStats.D_songs{Indices(j)}))];
        end
    else
        BirdStats.ConditionSpecific_DSong_INs{i} = [];
        BirdStats.ConditionSpecific_DSong_NumMotifs{i} = [];
        BirdStats.ConditionSpecific_DSong_FirstMotifDur{i} = [];
        
        BirdStats.ConditionSpecific_DSong_MeanINs{i} = [];
        BirdStats.ConditionSpecific_DSong_MeanNumMotifs{i} = [];
        BirdStats.ConditionSpecific_DSong_MeanFirstMotifDur{i} = [];

        for j = 1:length(Indices),
            BirdStats.ConditionSpecific_DSong_INs{i} = [BirdStats.ConditionSpecific_DSong_INs{i} BirdStats.BirdParameters(Indices(j)).NumINs(BirdStats.UN_songs{Indices(j)})];
            BirdStats.AllValidNumINs = [BirdStats.AllValidNumINs; [BirdStats.BirdParameters(Indices(j)).NumINs(BirdStats.UN_songs{Indices(j)})' ones(length(BirdStats.UN_songs{Indices(j)}),1)*BirdStats.ConditionIndex(Indices(j)) ones(length(BirdStats.UN_songs{Indices(j)}),1)*Indices(j) ones(length(BirdStats.UN_songs{Indices(j)}),1)*BirdStats.RecordingDayIndex(Indices(j)) ones(length(BirdStats.UN_songs{Indices(j)}),1)*BirdStats.BirdId]];
            
            BirdStats.ConditionSpecific_DSong_NumMotifs{i} = [BirdStats.ConditionSpecific_DSong_NumMotifs{i} BirdStats.BirdParameters(Indices(j)).NumMotifs(BirdStats.UN_songs{Indices(j)})];
            BirdStats.AllValidNumMotifs = [BirdStats.AllValidNumMotifs; [BirdStats.BirdParameters(Indices(j)).NumMotifs(BirdStats.UN_songs{Indices(j)})' ones(length(BirdStats.UN_songs{Indices(j)}),1)*BirdStats.ConditionIndex(Indices(j)) ones(length(BirdStats.UN_songs{Indices(j)}),1)*Indices(j) ones(length(BirdStats.UN_songs{Indices(j)}),1)*BirdStats.RecordingDayIndex(Indices(j)) ones(length(BirdStats.UN_songs{Indices(j)}),1)*BirdStats.BirdId]];
            
            BirdStats.ConditionSpecific_DSong_FirstMotifDur{i} = [BirdStats.ConditionSpecific_DSong_FirstMotifDur{i} BirdStats.BirdParameters(Indices(j)).FirstMotifDur(BirdStats.UN_songs{Indices(j)})];
            BirdStats.AllValidFirstMotifDur = [BirdStats.AllValidFirstMotifDur; [BirdStats.BirdParameters(Indices(j)).FirstMotifDur(BirdStats.UN_songs{Indices(j)})' ones(length(BirdStats.UN_songs{Indices(j)}),1)*BirdStats.ConditionIndex(Indices(j)) ones(length(BirdStats.UN_songs{Indices(j)}),1)*Indices(j) ones(length(BirdStats.UN_songs{Indices(j)}),1)*BirdStats.RecordingDayIndex(Indices(j)) ones(length(BirdStats.UN_songs{Indices(j)}),1)*BirdStats.BirdId]];
            
            BirdStats.ConditionSpecific_DSong_MeanINs{i} = [BirdStats.ConditionSpecific_DSong_MeanINs{i} mean(BirdStats.BirdParameters(Indices(j)).NumINs(BirdStats.UN_songs{Indices(j)}))];
            BirdStats.ConditionSpecific_DSong_MeanNumMotifs{i} = [BirdStats.ConditionSpecific_DSong_MeanNumMotifs{i} mean(BirdStats.BirdParameters(Indices(j)).NumMotifs(BirdStats.UN_songs{Indices(j)}))];
            BirdStats.ConditionSpecific_DSong_MeanFirstMotifDur{i} = [BirdStats.ConditionSpecific_DSong_MeanFirstMotifDur{i} nanmean(BirdStats.BirdParameters(Indices(j)).FirstMotifDur(BirdStats.UN_songs{Indices(j)}))];
            
        end
    end
end

% For each bird plot the # of INs, box plots for D songs across days and UN
% songs for undirected condition.
% And then plot the mean per condition as averages for each day
% And then plot the mean per condition pooled

close all
figure;
p = panel();
p.pack({1/2 1/2});
p(2).pack('h', {1/2 1/2});

p(1).select();
for i = 1:length(BirdStats.Date),
    if (strmatch('UN', BirdStats.Condition{i}))
        if (~isempty(BirdStats.UN_songs{i}))
            plot(i, BirdStats.BirdParameters(i).NumINs(BirdStats.UN_songs{i}), 'ko');
            hold on;
            BirdStats.MeanINNum(i) = mean(BirdStats.BirdParameters(i).NumINs(BirdStats.UN_songs{i}));
            errorbar(i-0.2, mean(BirdStats.BirdParameters(i).NumINs(BirdStats.UN_songs{i})), std(BirdStats.BirdParameters(i).NumINs(BirdStats.UN_songs{i}))/sqrt(length(BirdStats.UN_songs{i})), 'ks-', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
        else
            BirdStats.MeanINNum(i) = NaN;
        end
    else
        if (~isempty(BirdStats.D_songs{i}))
            plot(i, BirdStats.BirdParameters(i).NumINs(BirdStats.D_songs{i}), 'ko');
            hold on;
            errorbar(i-0.2, mean(BirdStats.BirdParameters(i).NumINs(BirdStats.D_songs{i})), std(BirdStats.BirdParameters(i).NumINs(BirdStats.D_songs{i}))/sqrt(length(BirdStats.D_songs{i})), 'ks-', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
            BirdStats.MeanINNum(i) = mean(BirdStats.BirdParameters(i).NumINs(BirdStats.D_songs{i}));
        else
            BirdStats.MeanINNum(i) = NaN;
        end
    end
    XLabelString{i} = [BirdStats.Date{i}, ' : ', BirdStats.Condition{i}];
end
axis tight;
Temp = axis;
Temp = [0.5 length(BirdStats.Date)+0.5 0 1.1*Temp(4)];
axis(Temp);
set(gca, 'XTick', 1:1:length(BirdStats.Date), 'XTickLabel', XLabelString, 'XTickLabelRotation', 30);
ylabel('# of INs');
title(BirdStats.BirdName);

p(2,1).select();
% Now plot IN number as averages of different sessions for each condition
for i = 1:length(BirdStats.UniqueConditions),
    plot(i+rand(length(BirdStats.ConditionSpecific_DSong_MeanINs{i}),1)/10, BirdStats.ConditionSpecific_DSong_MeanINs{i}, 'ko');
    hold on;
    errorbar(i-0.2, nanmean(BirdStats.ConditionSpecific_DSong_MeanINs{i}), nanstd(BirdStats.ConditionSpecific_DSong_MeanINs{i})/sqrt(length(find(~isnan(BirdStats.ConditionSpecific_DSong_MeanINs{i})))), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 8); 
end
set(gca, 'XTick', 1:1:length(BirdStats.UniqueConditions), 'XTickLabel', BirdStats.UniqueConditions);
axis tight;
Temp = axis;
Temp = [0.5 length(BirdStats.UniqueConditions)+0.5 0 1.1*Temp(4)];
axis(Temp);
ylabel('# of INs');
xlabel('Distance from the female');

p(2,2).select();
% Now plot IN number as averages of all numbers pooled together
for i = 1:length(BirdStats.UniqueConditions),
    plot(i+rand(length(BirdStats.ConditionSpecific_DSong_INs{i}),1)/10, BirdStats.ConditionSpecific_DSong_INs{i}, 'ko');
    hold on;
    errorbar(i-0.2, nanmean(BirdStats.ConditionSpecific_DSong_INs{i}), nanstd(BirdStats.ConditionSpecific_DSong_INs{i})/sqrt(length(find(~isnan(BirdStats.ConditionSpecific_DSong_INs{i})))), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 8); 
end
set(gca, 'XTick', 1:1:length(BirdStats.UniqueConditions), 'XTickLabel', BirdStats.UniqueConditions);
axis tight;
Temp2 = axis;
Temp2 = [0.5 length(BirdStats.UniqueConditions)+0.5 0 max(1.1*Temp2(4), Temp(4))];
axis(Temp2);
xlabel('Distance from the female');

p.margintop = 10;
p.fontsize = 14;
p.marginleft = 20;
p.de.margin = 25;
p(2,1).select();
axis(Temp2);
set(gcf, 'Color', 'w');
set(gcf, 'Position', [571 264 1150 550]);
OutputDir = '/home/raghav/StudentRelated/Harini/Manuscript/JEBManuscript/Figures/';
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(OutputDir, [BirdStats.BirdName, '.NumINs.png']), '-dpng');

figure;
p = panel();
p.pack({1/2 1/2});
p(2).pack('h', {1/2 1/2});

p(1).select();
for i = 1:length(BirdStats.Date),
    if (strmatch('UN', BirdStats.Condition{i}))
        if (~isempty(BirdStats.UN_songs{i}))
            plot(i, BirdStats.BirdParameters(i).NumMotifs(BirdStats.UN_songs{i}), 'ko');
            hold on;
            BirdStats.MeanMotifNum(i) = mean(BirdStats.BirdParameters(i).NumMotifs(BirdStats.UN_songs{i}));
            errorbar(i-0.2, mean(BirdStats.BirdParameters(i).NumMotifs(BirdStats.UN_songs{i})), std(BirdStats.BirdParameters(i).NumMotifs(BirdStats.UN_songs{i}))/sqrt(length(BirdStats.UN_songs{i})), 'ks-', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
        else
            BirdStats.MeanMotifNum(i) = NaN;
        end
    else
        if (~isempty(BirdStats.D_songs{i}))
            plot(i, BirdStats.BirdParameters(i).NumMotifs(BirdStats.D_songs{i}), 'ko');
            hold on;
            errorbar(i-0.2, mean(BirdStats.BirdParameters(i).NumMotifs(BirdStats.D_songs{i})), std(BirdStats.BirdParameters(i).NumMotifs(BirdStats.D_songs{i}))/sqrt(length(BirdStats.D_songs{i})), 'ks-', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
            BirdStats.MeanMotifNum(i) = mean(BirdStats.BirdParameters(i).NumMotifs(BirdStats.D_songs{i}));
        else
            BirdStats.MeanMotifNum(i) = NaN;
        end
    end
    XLabelString{i} = [BirdStats.Date{i}, ' : ', BirdStats.Condition{i}];
end
axis tight;
Temp = axis;
Temp = [0.5 length(BirdStats.Date)+0.5 0 1.1*Temp(4)];
axis(Temp);
set(gca, 'XTick', 1:1:length(BirdStats.Date), 'XTickLabel', XLabelString, 'XTickLabelRotation', 30);
ylabel('# of motifs/bout');
title(BirdStats.BirdName);

p(2,1).select();
% Now plot IN number as averages of different sessions for each condition
for i = 1:length(BirdStats.UniqueConditions),
    plot(i+rand(length(BirdStats.ConditionSpecific_DSong_MeanNumMotifs{i}),1)/10, BirdStats.ConditionSpecific_DSong_MeanNumMotifs{i}, 'ko');
    hold on;
    errorbar(i-0.2, nanmean(BirdStats.ConditionSpecific_DSong_MeanNumMotifs{i}), nanstd(BirdStats.ConditionSpecific_DSong_MeanNumMotifs{i})/sqrt(length(find(~isnan(BirdStats.ConditionSpecific_DSong_MeanNumMotifs{i})))), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 8); 
end
set(gca, 'XTick', 1:1:length(BirdStats.UniqueConditions), 'XTickLabel', BirdStats.UniqueConditions);
axis tight;
Temp = axis;
Temp = [0.5 length(BirdStats.UniqueConditions)+0.5 0 1.1*Temp(4)];
axis(Temp);
ylabel('# of motifs/bout');
xlabel('Distance from the female');

p(2,2).select();
% Now plot IN number as averages of all numbers pooled together
for i = 1:length(BirdStats.UniqueConditions),
    plot(i+rand(length(BirdStats.ConditionSpecific_DSong_NumMotifs{i}),1)/10, BirdStats.ConditionSpecific_DSong_NumMotifs{i}, 'ko');
    hold on;
    errorbar(i-0.2, nanmean(BirdStats.ConditionSpecific_DSong_NumMotifs{i}), nanstd(BirdStats.ConditionSpecific_DSong_NumMotifs{i})/sqrt(length(find(~isnan(BirdStats.ConditionSpecific_DSong_NumMotifs{i})))), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 8); 
end
set(gca, 'XTick', 1:1:length(BirdStats.UniqueConditions), 'XTickLabel', BirdStats.UniqueConditions);
axis tight;
Temp2 = axis;
Temp2 = [0.5 length(BirdStats.UniqueConditions)+0.5 0 max(1.1*Temp2(4), Temp(4))];
axis(Temp2);
xlabel('Distance from the female');

p.margintop = 10;
p.fontsize = 14;
p.marginleft = 20;
p.de.margin = 25;
p(2,1).select();
axis(Temp2);
set(gcf, 'Color', 'w');
set(gcf, 'Position', [571 264 1150 550]);
OutputDir = '/home/raghav/StudentRelated/Harini/Manuscript/JEBManuscript/Figures/';
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(OutputDir, [BirdStats.BirdName, '.NumMotifs.png']), '-dpng');

figure;
p = panel();
p.pack({1/2 1/2});
p(2).pack('h', {1/2 1/2});

p(1).select();
for i = 1:length(BirdStats.Date),
    if (strmatch('UN', BirdStats.Condition{i}))
        if (~isempty(BirdStats.UN_songs{i}))
            plot(i, BirdStats.BirdParameters(i).FirstMotifDur(BirdStats.UN_songs{i}), 'ko');
            hold on;
            BirdStats.MeanFirstMotifDur(i) = nanmean(BirdStats.BirdParameters(i).FirstMotifDur(BirdStats.UN_songs{i}));
            errorbar(i-0.2, nanmean(BirdStats.BirdParameters(i).FirstMotifDur(BirdStats.UN_songs{i})), nanstd(BirdStats.BirdParameters(i).FirstMotifDur(BirdStats.UN_songs{i}))/sqrt(length(find(~isnan(BirdStats.BirdParameters(i).FirstMotifDur(BirdStats.UN_songs{i}))))), 'ks-', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
        else
            BirdStats.MeanFirstMotifDur(i) = NaN;
        end
    else
        if (~isempty(BirdStats.D_songs{i}))
            plot(i, BirdStats.BirdParameters(i).FirstMotifDur(BirdStats.D_songs{i}), 'ko');
            hold on;
            errorbar(i-0.2, nanmean(BirdStats.BirdParameters(i).FirstMotifDur(BirdStats.D_songs{i})), nanstd(BirdStats.BirdParameters(i).FirstMotifDur(BirdStats.D_songs{i}))/sqrt(length(find(~isnan(BirdStats.BirdParameters(i).FirstMotifDur(BirdStats.D_songs{i}))))), 'ks-', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
            BirdStats.MeanFirstMotifDur(i) = nanmean(BirdStats.BirdParameters(i).FirstMotifDur(BirdStats.D_songs{i}));
        else
            BirdStats.MeanINNum(i) = NaN;
        end
    end
    XLabelString{i} = [BirdStats.Date{i}, ' : ', BirdStats.Condition{i}];
end
axis tight;
Temp = axis;
Temp = [0.5 length(BirdStats.Date)+0.5 0.95*Temp(3) 1.05*Temp(4)];
axis(Temp);
set(gca, 'XTick', 1:1:length(BirdStats.Date), 'XTickLabel', XLabelString, 'XTickLabelRotation', 30);
ylabel('First motif duration (msec)');
title(BirdStats.BirdName);

p(2,1).select();
% Now plot IN number as averages of different sessions for each condition
for i = 1:length(BirdStats.UniqueConditions),
    plot(i+rand(length(BirdStats.ConditionSpecific_DSong_MeanFirstMotifDur{i}),1)/10, BirdStats.ConditionSpecific_DSong_MeanFirstMotifDur{i}, 'ko');
    hold on;
    errorbar(i-0.2, nanmean(BirdStats.ConditionSpecific_DSong_MeanFirstMotifDur{i}), nanstd(BirdStats.ConditionSpecific_DSong_MeanFirstMotifDur{i})/sqrt(length(find(~isnan(BirdStats.ConditionSpecific_DSong_MeanFirstMotifDur{i})))), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 8); 
end
set(gca, 'XTick', 1:1:length(BirdStats.UniqueConditions), 'XTickLabel', BirdStats.UniqueConditions);
axis tight;
Temp = axis;
Temp = [0.5 length(BirdStats.UniqueConditions)+0.5 0.95*Temp(3) 1.05*Temp(4)];
axis(Temp);
ylabel('First motif duration (msec)');
xlabel('Distance from the female');

p(2,2).select();
% Now plot IN number as averages of all numbers pooled together
for i = 1:length(BirdStats.UniqueConditions),
    plot(i+rand(length(BirdStats.ConditionSpecific_DSong_FirstMotifDur{i}),1)/10, BirdStats.ConditionSpecific_DSong_FirstMotifDur{i}, 'ko');
    hold on;
    errorbar(i-0.2, nanmean(BirdStats.ConditionSpecific_DSong_FirstMotifDur{i}), nanstd(BirdStats.ConditionSpecific_DSong_INs{i})/sqrt(length(find(~isnan(BirdStats.ConditionSpecific_DSong_FirstMotifDur{i})))), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 8); 
end
set(gca, 'XTick', 1:1:length(BirdStats.UniqueConditions), 'XTickLabel', BirdStats.UniqueConditions);
axis tight;
Temp2 = axis;
Temp2 = [0.5 length(BirdStats.UniqueConditions)+0.5 min(0.95*Temp2(3), Temp(3)) max(1.05*Temp2(4), Temp(4))];
axis(Temp2);
xlabel('Distance from the female');

p.de.margin = 25;
p.margintop = 10;
p(2,1).select();
axis(Temp2);
set(gcf, 'Color', 'w');
set(gcf, 'Position', [571 264 1150 750]);
p.fontsize = 14;
p.marginleft = 20;
OutputDir = '/home/raghav/StudentRelated/Harini/Manuscript/JEBManuscript/Figures/';
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(OutputDir, [BirdStats.BirdName, '.FirstMotifDur.png']), '-dpng');

disp('Finished');