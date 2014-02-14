function [] = PlotAllSyllSAPFeats(SyllFeats, MotifLengths, Syllables, PlotMarkers, BirdName)

FeatureCols{1} = 'Duration';
FeatureCols{2} = 'Amplitude';
FeatureCols{3} = 'Entropy';
FeatureCols{4} = 'Mean Frequency';
FeatureCols{5} = 'Amplitude Modulation';
FeatureCols{6} = 'Pitch Goodness';
FeatureCols{7} = 'Frequency Modulation';

figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [400 100 900 600]);

for i = 1:length(Syllables),
    for j = 1:length(FeatureCols),
        subplot(2,4,j);
        title(FeatureCols{j});
        hold on;
        Index1 = -100;
        Index2 = -100;
        for k = 1:length(SyllFeats{1}),
            if (SyllFeats{1}{k}.Label == Syllables(i))
                Index1 = k;
                break;
            end
        end
        for k = 1:length(SyllFeats{2}),
            if (SyllFeats{2}{k}.Label == Syllables(i))
                Index2 = k;
                break;
            end
        end
        
        if ((Index1 > 0) && (Index2 > 0))
            plot(mean(SyllFeats{1}{Index1}.Feats(:,j)), mean(SyllFeats{2}{Index2}.Feats(:,j)), PlotMarkers(i,:));
            plot([mean(SyllFeats{1}{Index1}.Feats(:,j)) mean(SyllFeats{1}{Index1}.Feats(:,j))], [(mean(SyllFeats{2}{Index2}.Feats(:,j)) + std(SyllFeats{2}{Index2}.Feats(:,j))) (mean(SyllFeats{2}{Index2}.Feats(:,j)) - std(SyllFeats{2}{Index2}.Feats(:,j)))], PlotMarkers(i,1));
            plot([(mean(SyllFeats{1}{Index1}.Feats(:,j)) + std(SyllFeats{1}{Index1}.Feats(:,j))) (mean(SyllFeats{1}{Index1}.Feats(:,j)) - std(SyllFeats{1}{Index1}.Feats(:,j)))], [mean(SyllFeats{2}{Index2}.Feats(:,j)) mean(SyllFeats{2}{Index2}.Feats(:,j))], PlotMarkers(i,1));
        end
    end
end

for i = 1:length(FeatureCols),
    subplot(2,4,i);
    axis tight;
    temp = axis;
    axis([min([temp(1) temp(3)]) max([temp(2) temp(4)]) min([temp(1) temp(3)]) max([temp(2) temp(4)])]);
    plot([min([temp(1) temp(3)]) max([temp(2) temp(4)])], [min([temp(1) temp(3)]) max([temp(2) temp(4)])], 'k--');
    if (i > 3)
        xlabel('Pre');
    end
    if (mod(i,4) == 1)
        ylabel('Post 1');
    end
    set(gca, 'FontSize', 12);
end

subplot(2,4,8);
PreBar = bar(1, mean(MotifLengths{1}));
set(PreBar, 'FaceColor', 'none', 'EdgeColor', 'k');
hold on;
plot([1 1], [(mean(MotifLengths{1}) + std(MotifLengths{1})) (mean(MotifLengths{1}) - std(MotifLengths{1}))], 'k');

Post1Bar = bar(2, mean(MotifLengths{2}));
set(Post1Bar, 'FaceColor', 'none', 'EdgeColor', 'k');
hold on;
plot([2 2], [(mean(MotifLengths{2}) + std(MotifLengths{2})) (mean(MotifLengths{2}) - std(MotifLengths{2}))], 'k');

axis tight;
Temp = axis;
Temp(3) = 0;
Temp(1:2) = [0.5 2.5];
Temp(4) = 1.1*Temp(4);
axis(Temp);
set(gca, 'Box', 'off');
title('Motif Durations');
set(gca, 'XTick', [1 2], 'XTickLabel', [{'Pre'} {'Post 1'}]);
ylabel('Motif length (sec)');

set(gcf, 'PaperPositionMode', 'auto');
print(gcf, '-opengl', '-dtiff', '-r0', [BirdName, '_Dir_SAPFeats_Pre_Post1.tif']);

        
% Now everything is plotted normalized to the mean of Pre
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [400 100 900 600]);

for i = 1:length(Syllables),
    for j = 1:length(FeatureCols),
        subplot(2,4,j);
        title(FeatureCols{j});
        hold on;
        Index1 = -100;
        Index2 = -100;
        for k = 1:length(SyllFeats{1}),
            if (SyllFeats{1}{k}.Label == Syllables(i))
                Index1 = k;
                break;
            end
        end
        for k = 1:length(SyllFeats{2}),
            if (SyllFeats{2}{k}.Label == Syllables(i))
                Index2 = k;
                break;
            end
        end
        
        if ((Index1 > 0) && (Index2 > 0))
            % normalizing it to the mean of Pre for that syllable
            SyllFeats{2}{Index2}.Feats(:,j) = SyllFeats{2}{Index2}.Feats(:,j)/mean(SyllFeats{1}{Index1}.Feats(:,j));            
            SyllFeats{1}{Index1}.Feats(:,j) = SyllFeats{1}{Index1}.Feats(:,j)/mean(SyllFeats{1}{Index1}.Feats(:,j));
                        
            plot(i, mean(SyllFeats{2}{Index2}.Feats(:,j)), PlotMarkers(i,:));
            plot([i i], [(mean(SyllFeats{2}{Index2}.Feats(:,j)) + std(SyllFeats{2}{Index2}.Feats(:,j))) (mean(SyllFeats{2}{Index2}.Feats(:,j)) - std(SyllFeats{2}{Index2}.Feats(:,j)))], PlotMarkers(i,1));
            %plot([(mean(SyllFeats{1}{Index1}.Feats(:,j)) + std(SyllFeats{1}{Index1}.Feats(:,j))) (mean(SyllFeats{1}{Index1}.Feats(:,j)) - std(SyllFeats{1}{Index1}.Feats(:,j)))], [mean(SyllFeats{2}{Index2}.Feats(:,j)) mean(SyllFeats{2}{Index2}.Feats(:,j))], PlotMarkers(i,1));
        end
    end
end

for i = 1:length(Syllables),
    XLabels{i} = Syllables(i);
end
for i = 1:length(FeatureCols),
    subplot(2,4,i);
    axis tight;
    temp = axis;
    axis([0 length(Syllables)+1 0.95*temp(3) 1.05*temp(4)]);
    plot([0 length(Syllables)+1], [1 1], 'k--');
    set(gca, 'XTick', [1:1:length(Syllables)], 'XTickLabel', XLabels);
    if (mod(i,4) == 1)
        ylabel('Post1/Pre');
    end
    set(gca, 'FontSize', 12);
end

subplot(2,4,8);
% normalizing motif lengths to mean of Pre
MotifLengths{2} = MotifLengths{2}/mean(MotifLengths{1});
MotifLengths{1} = MotifLengths{1}/mean(MotifLengths{1});


PreBar = bar(1, mean(MotifLengths{1}));
set(PreBar, 'FaceColor', 'none', 'EdgeColor', 'k');
hold on;
plot([1 1], [(mean(MotifLengths{1}) + std(MotifLengths{1})) (mean(MotifLengths{1}) - std(MotifLengths{1}))], 'k');

Post1Bar = bar(2, mean(MotifLengths{2}));
set(Post1Bar, 'FaceColor', 'none', 'EdgeColor', 'k');
hold on;
plot([2 2], [(mean(MotifLengths{2}) + std(MotifLengths{2})) (mean(MotifLengths{2}) - std(MotifLengths{2}))], 'k');

axis tight;
Temp = axis;
Temp(3) = 0;
Temp(1:2) = [0.5 2.5];
Temp(4) = 1.1*Temp(4);
axis(Temp);
set(gca, 'Box', 'off');
title('Motif Durations');
set(gca, 'XTick', [1 2], 'XTickLabel', [{'Pre'} {'Post 1'}]);
ylabel('Motif length (sec)');

set(gcf, 'PaperPositionMode', 'auto');
print(gcf, '-opengl', '-dtiff', '-r0', [BirdName, '_Dir_Normalized_SAPFeats_Pre_Post1.tif']);
