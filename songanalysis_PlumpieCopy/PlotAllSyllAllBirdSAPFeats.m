function [] = PlotAllSyllAllBirdSAPFeats(SyllFeats, Syllables, RemainingHVC)

FeatureCols{1} = 'Duration';
FeatureCols{2} = 'Amplitude';
FeatureCols{3} = 'Entropy';
FeatureCols{4} = 'Mean Frequency';
FeatureCols{5} = 'Amplitude Modulation';
FeatureCols{6} = 'Pitch Goodness';
FeatureCols{7} = 'Frequency Modulation';

Colours = 'rgbcm';
Symbols = 'od^';

figure(1);
set(gcf, 'Color', 'w');
set(gcf, 'Position', [400 100 900 600]);

figure(2);
set(gcf, 'Color', 'w');
set(gcf, 'Position', [400 100 900 600]);

AllBirdSyllables = [];
for BirdIndex = 1:length(SyllFeats),
    BirdSyllables = [];
    for i = 1:length(Syllables{BirdIndex}),
        for j = 1:4,
            figure(1);
            subplot(2,2,j);
            title(FeatureCols{j});
            hold on;
            Index1 = -100;
            Index2 = -100;
            for k = 1:length(SyllFeats{BirdIndex}{1}),
                if (SyllFeats{BirdIndex}{1}{k}.Label == Syllables{BirdIndex}(i))
                    Index1 = k;
                    break;
                end
            end
            for k = 1:length(SyllFeats{BirdIndex}{2}),
                if (SyllFeats{BirdIndex}{2}{k}.Label == Syllables{BirdIndex}(i))
                    Index2 = k;
                    break;
                end
            end

            if ((Index1 > 0) && (Index2 > 0))
                plot(mean(SyllFeats{BirdIndex}{1}{Index1}.Feats(:,j)), mean(SyllFeats{BirdIndex}{2}{Index2}.Feats(:,j)), [Colours(mod(BirdIndex-1, length(Colours)) + 1), Symbols(ceil(BirdIndex/length(Colours)))], 'MarkerSize', 4);
                figure(2);
                subplot(2,2,j);
                hold on;
                plot(RemainingHVC(BirdIndex), (mean(SyllFeats{BirdIndex}{2}{Index2}.Feats(:,j)) - mean(SyllFeats{BirdIndex}{1}{Index1}.Feats(:,j)))/mean(SyllFeats{BirdIndex}{1}{Index1}.Feats(:,j)) * 100,  [Colours(mod(BirdIndex-1, length(Colours)) + 1), Symbols(ceil(BirdIndex/length(Colours)))]);
                BirdSyllables(i,j) = (mean(SyllFeats{BirdIndex}{2}{Index2}.Feats(:,j)) - mean(SyllFeats{BirdIndex}{1}{Index1}.Feats(:,j)))/mean(SyllFeats{BirdIndex}{1}{Index1}.Feats(:,j)) * 100;
%                plot([mean(SyllFeats{BirdIndex}{1}{Index1}.Feats(:,j)) mean(SyllFeats{BirdIndex}{1}{Index1}.Feats(:,j))], [(mean(SyllFeats{BirdIndex}{2}{Index2}.Feats(:,j)) + std(SyllFeats{BirdIndex}{2}{Index2}.Feats(:,j))) (mean(SyllFeats{BirdIndex}{2}{Index2}.Feats(:,j)) - std(SyllFeats{BirdIndex}{2}{Index2}.Feats(:,j)))], 'k');
%                plot([(mean(SyllFeats{BirdIndex}{1}{Index1}.Feats(:,j)) + std(SyllFeats{BirdIndex}{1}{Index1}.Feats(:,j))) (mean(SyllFeats{BirdIndex}{1}{Index1}.Feats(:,j)) - std(SyllFeats{BirdIndex}{1}{Index1}.Feats(:,j)))], [mean(SyllFeats{BirdIndex}{2}{Index2}.Feats(:,j)) mean(SyllFeats{BirdIndex}{2}{Index2}.Feats(:,j))], 'k');
            end
        end
    end
    for j = 1:4,
        figure(2);
        subplot(2,2,j);
        errorbar(RemainingHVC(BirdIndex), mean(BirdSyllables(:,j)), std(BirdSyllables(:,j)), ['ks'], 'MarkerSize', 6);
    end
    AllBirdSyllables = [AllBirdSyllables; [ones(size(BirdSyllables(:,1)))*RemainingHVC(BirdIndex) BirdSyllables]];
end

for i = 1:4,
    figure(1);
    subplot(2,2,i);
    axis tight;
    temp = axis;
    axis([min([temp(1) temp(3)]) max([temp(2) temp(4)]) min([temp(1) temp(3)]) max([temp(2) temp(4)])]);
    plot([min([temp(1) temp(3)]) max([temp(2) temp(4)])], [min([temp(1) temp(3)]) max([temp(2) temp(4)])], 'k--');
    if (i > 2)
        xlabel('Pre');
    end
    if (mod(i,2) == 1)
        ylabel('Post 1');
    end
    set(gca, 'FontSize', 12);
    figure(2);
    subplot(2,2,i);
    ylabel(['% change in ', FeatureCols{i}], 'FontSize', 14);
    set(gca, 'Fontsize', 12);
    if (i > 2)
        xlabel('% HVC remaining', 'FontSize', 14);
    end
    axis tight;
    temp = axis;
    temp = [0.98*temp(1) 1.02*temp(2) 0.98*temp(3) 1.02*temp(4)];
    axis(temp);
    plot([temp(1) temp(2)], [0 0], 'k--');
end

figure(2);
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, '-opengl', '-dtiff', '-r0', 'AllBirds_Dir_SAPFeats_Pre_Post1.tif');

       