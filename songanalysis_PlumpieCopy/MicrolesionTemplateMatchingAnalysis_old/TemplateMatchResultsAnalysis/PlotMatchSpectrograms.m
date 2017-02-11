function [] = PlotMatchSpectrograms(handles, Syllable, varargin)


% File to plot the spectrograms of template match values for each of the
% days

if (nargin > 2)
    InputData = varargin{1};
end

Colours = 'rgbcmky';

PlotFigure = figure;
set(PlotFigure, 'Position', [427 170 675 525]);

SyllIndex = intersect(find(cellfun(@length, strfind([handles(1).Labels], Syllable))), find(cellfun(@length, handles(1).Labels) == length(Syllable)));
for i = 1:length(handles),
    Legend{i} = ['Day #', num2str(i)];
 
    % First analyse the output data from directed songs
    
    if (~isempty(handles(i).DirPeaks{SyllIndex}))
        MotifPeaks = handles(i).DirPeaks{SyllIndex}(:,1);
        DirMotifStats = handles(i).DirPeaks{SyllIndex};
    
        figure(PlotFigure);
        subplot(2, 1, 1);
        plot(0:0.1:max(MotifPeaks), 100*histc(MotifPeaks, (0:0.1:max(MotifPeaks)))/sum(histc(MotifPeaks, (0:0.1:max(MotifPeaks)))), Colours(i), 'LineWidth', 2);
        hold on;
    end
    
    % Next analyse the output data from undirected songs
    MotifPeaks = handles(i).UnDirPeaks{SyllIndex}(:,1);
    UnDirMotifStats = handles(i).UnDirPeaks{SyllIndex};
    figure(PlotFigure);
    subplot(2, 1, 2);
    plot(0:0.1:max(MotifPeaks), 100*histc(MotifPeaks, (0:0.1:max(MotifPeaks)))/sum(histc(MotifPeaks, (0:0.1:max(MotifPeaks)))), Colours(i), 'LineWidth', 2);
    hold on;
end

figure(PlotFigure);
set(gcf, 'Color', 'w');
subplot(2,1,1);
axis tight;
TempAxis(1,:) = axis;
    
subplot(2,1,2);
TempAxis(2,:) = axis;
axis tight;

subplot(2,1,1);
set(gca, 'Box', 'off');
axis([0 max(TempAxis(:,2))*1.05 0 0.5]);
if (length(Syllable) > 1)
    title(['Directed motif ', Syllable], 'FontSize', 12, 'FontWeight', 'bold');
else
    title(['Directed syllable ', Syllable], 'FontSize', 12, 'FontWeight', 'bold');
end
ylabel('%', 'FontSize', 12, 'FontWeight', 'bold');
legend(Legend);
set(gca, 'FontSize', 12);

subplot(2,1,2);
set(gca, 'Box', 'off');
axis([0 max(TempAxis(:,2))*1.05 0 0.5]);
if (length(Syllable) > 1)
    title(['UnDirected motif ', Syllable], 'FontSize', 12, 'FontWeight', 'bold');
else
    title(['UnDirected syllable ', Syllable], 'FontSize', 12, 'FontWeight', 'bold');
end
xlabel('Template match value', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('%', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12);

% Now for all the random comparisons
figure(PlotFigure);
Criterion = mean(handles(1).RandomSongComparisonUnDirMaxPeaks{SyllIndex}(:,1)) + 3*std(handles(1).RandomSongComparisonUnDirMaxPeaks{SyllIndex}(:,1));
subplot(2, 1, 1);
plot([Criterion Criterion], [0 50], [Colours(1), ':'], 'LineWidth', 2);
subplot(2, 1, 2);
plot([Criterion Criterion], [0 50], [Colours(1), ':'], 'LineWidth', 2);
hold on;

figure(PlotFigure);
Criterion = mean(handles(1).RandomSoundComparisonUnDirMaxPeaks{SyllIndex}(:,1)) + 3*std(handles(1).RandomSoundComparisonUnDirMaxPeaks{SyllIndex}(:,1));
subplot(2, 1, 1);
plot([Criterion Criterion], [0 50], [Colours(1), '--'], 'LineWidth', 2);
subplot(2, 1, 2);
plot([Criterion Criterion], [0 50], [Colours(1), '--'], 'LineWidth', 2);
hold on;

figure(PlotFigure);
Criterion = mean(handles(1).ShuffledSongComparisonsUnDirMaxPeaks{SyllIndex}(:,1)) + 3*std(handles(1).ShuffledSongComparisonsUnDirMaxPeaks{SyllIndex}(:,1));
subplot(2, 1, 1);
plot([Criterion Criterion], [0 50], [Colours(1), '-.'], 'LineWidth', 2);
subplot(2, 1, 2);
plot([Criterion Criterion], [0 50], [Colours(1), '-.'], 'LineWidth', 2);
hold on;

if (SyllIndex > 1)
    figure(PlotFigure);
    for i = 1:length(handles(1).UnDirSyllMaxPeaks),
        if ((i+1) ~= SyllIndex)
            NotNaNIndices = ~isnan(handles(1).UnDirSyllMaxOtherPeaks{SyllIndex-1}{i}(:,1));
            Criterion = mean(handles(1).UnDirSyllMaxOtherPeaks{SyllIndex-1}{i}(NotNaNIndices,1)) + 3*std(handles(1).UnDirSyllMaxOtherPeaks{SyllIndex-1}{i}(NotNaNIndices,1));
            subplot(2, 1, 1);
            plot([Criterion Criterion], [0 50], [Colours(2), '-.'], 'LineWidth', 2);
            subplot(2, 1, 2);
            plot([Criterion Criterion], [0 50], [Colours(2), '-.'], 'LineWidth', 2);
            hold on;
            text(Criterion, 0.4, handles(1).Labels{i+1}, 'FontSize', 12);
        end
    end
end

Flag = 1;

if (~isempty(strfind(InputData.MotifTemplate.MotifTemplate(1).Label, Syllable)) && (length(Syllable) == length(InputData.MotifTemplate.MotifTemplate(1).Label)))
    ZeroIndex = find(([InputData.MotifTemplate.MotifTemplate.TimeStretch] == 0) & ([InputData.MotifTemplate.MotifTemplate.FreqStretch] == 0));
    TemplateLen = size(InputData.MotifTemplate.MotifTemplate(ZeroIndex).MotifTemplate, 2)/250;
else
    for i = 1:length(InputData.SyllTemplates.SyllableTemplates),
        if (~isempty(find(InputData.SyllTemplates.SyllableTemplates{i}{1}.MotifTemplate(1).Label == Syllable)))
            ZeroIndex = find(([InputData.SyllTemplates.SyllableTemplates{i}{1}.MotifTemplate.TimeStretch] == 0) & ([InputData.SyllTemplates.SyllableTemplates{i}{1}.MotifTemplate.FreqStretch] == 0));
            TemplateLen = size(InputData.SyllTemplates.SyllableTemplates{i}{1}.MotifTemplate(ZeroIndex).MotifTemplate, 2)/250;
        end
    end
end

NoToPlot = 4;

while (Flag)
    figure(PlotFigure);
    annotation('textbox', [0.05 0.97 0.9 0.02], 'String', 'Choose a template match value to see representative spectrograms; type q to quit', 'EdgeColor', 'w');;
    [x, y, button] = ginput(1);
    if (button == 113)
        break;
    end
    figure;
    for i = 1:InputData.NoofDaysToAnalyse,
        if (~isempty(handles(i).DirPeaks{SyllIndex}))
            [Temp, TempSortedIndices] = sort(handles(i).DirPeaks{SyllIndex}(:,1));
            Matches = find(Temp >= x, 1, 'first');
            if ((isempty(Matches)) || ((Matches + (NoToPlot - 2)) > length(Temp)))
                Matches = TempSortedIndices(end-(NoToPlot-1):end);
            else
                Matches = TempSortedIndices((Matches - 1):(Matches + (NoToPlot - 2)));
            end
                
            for j = 1:length(Matches),
                subplot(8,InputData.NoofDaysToAnalyse,(InputData.NoofDaysToAnalyse)*(j-1) + i);
                PlotSpectrogramInAxis([InputData.DataDirectories{i}, '/'], handles(i).DirSongFiles{handles(i).DirPeaks{SyllIndex}(Matches(j),3)}, InputData.FileType, gca, [(handles(i).DirPeaks{SyllIndex}(Matches(j),2) - 0.2) (handles(i).DirPeaks{SyllIndex}(Matches(j),2) + TemplateLen + 0.2)]);
                hold on;
                plot([handles(i).DirPeaks{SyllIndex}(Matches(j),2) handles(i).DirPeaks{SyllIndex}(Matches(j),2)], [300 8000], 'k', 'LineWidth', 2);
                xlabel('');
                ylabel('');
                set(gca, 'YTickLabel', []);
                set(gca, 'XTick', [0 TemplateLen]);
                title(handles(i).DirPeaks{SyllIndex}(Matches(j),1));
            end
        end
        
        [Temp, TempSortedIndices] = sort(handles(i).UnDirPeaks{SyllIndex}(:,1));
        Matches = find(Temp >= x, 1, 'first');
        if ((isempty(Matches)) || ((Matches + (NoToPlot - 2)) > length(Temp)))
            Matches = TempSortedIndices(end-(NoToPlot-1):end);
        else
            Matches = TempSortedIndices((Matches - 1):(Matches + (NoToPlot - 2)));
        end
        for j = 1:length(Matches),
            subplot(8,InputData.NoofDaysToAnalyse,(InputData.NoofDaysToAnalyse)*(j-1+4) + i);
            PlotSpectrogramInAxis([InputData.DataDirectories{i}, '/'], handles(i).UnDirSongFiles{handles(i).UnDirPeaks{SyllIndex}(Matches(j),3)}, InputData.FileType, gca, [(handles(i).UnDirPeaks{SyllIndex}(Matches(j),2) - 0.2) (handles(i).UnDirPeaks{SyllIndex}(Matches(j),2) + TemplateLen + 0.2)]);
            plot([handles(i).UnDirPeaks{SyllIndex}(Matches(j),2) handles(i).UnDirPeaks{SyllIndex}(Matches(j),2)], [300 8000], 'k', 'LineWidth', 2);
            xlabel('');
            ylabel('');
            set(gca, 'YTickLabel', []);
            set(gca, 'XTick', [0 TemplateLen]);
            title(handles(i).UnDirPeaks{SyllIndex}(Matches(j),1));
        end
    end
    disp('Chose');
end
disp('Finished analysis');