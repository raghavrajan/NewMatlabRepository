function [Output] = PlotGoodMatchDetails(handles, Syllable, varargin)

% File to plot the spectrograms of template match values for each of the
% days

if (nargin > 2)
    InputData = varargin{1};
end

Colours = 'rgbcmky';

SyllIndex = intersect(find(cellfun(@length, strfind([handles(1).Labels], Syllable))), find(cellfun(@length, handles(1).Labels) == length(Syllable)));

% Now for all the random comparisons
Criterion(1) = mean(handles(1).RandomSongComparisonUnDirMaxPeaks{SyllIndex}(:,1)) + 3*std(handles(1).RandomSongComparisonUnDirMaxPeaks{SyllIndex}(:,1));
Criterion(2) = mean(handles(1).RandomSoundComparisonUnDirMaxPeaks{SyllIndex}(:,1)) + 3*std(handles(1).RandomSoundComparisonUnDirMaxPeaks{SyllIndex}(:,1));
Criterion(3) = mean(handles(1).ShuffledSongComparisonsUnDirMaxPeaks{SyllIndex}(:,1)) + 3*std(handles(1).ShuffledSongComparisonsUnDirMaxPeaks{SyllIndex}(:,1));

% Criterion = max(Criterion);
Criterion = Criterion(3);

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

PlotFigure = figure;
set(PlotFigure, 'Position', [427 170 675 525]);

for i = 1:length(handles),
    Legend{(i-1)*2 + 1} = ['Day #', num2str(i), ':D'];
 
    % First analyse the output data from directed songs
    
    if (~isempty(handles(i).DirPeaks{SyllIndex}))
        MotifPeaks = handles(i).DirPeaks{SyllIndex}(:,1);
        GoodPeakIndices = find(MotifPeaks >= Criterion);
    
        if (~isempty(GoodPeakIndices))
            figure(PlotFigure);
            subplot(3, 1, 1);
            errorbar(i-0.1, mean(MotifPeaks(GoodPeakIndices)), std(MotifPeaks(GoodPeakIndices)), 'rs');
            hold on;
            
            figure(PlotFigure);
            subplot(3, 1, 2);
            GoodPeaks = sort(MotifPeaks(GoodPeakIndices));
            plot(i-0.1, median(GoodPeaks), 'rs');
            hold on;
            plot([i-0.1 i-0.1], [GoodPeaks(ceil(0.25*length(GoodPeaks))) GoodPeaks(ceil(0.75*length(GoodPeaks)))], 'r');
        end
        
        figure(PlotFigure);
        subplot(3, 1, 3);
        TotalFileTime = sum(handles(i).DirFileLens);
        plot(i-0.1, length(GoodPeakIndices)*TemplateLen/TotalFileTime, 'rs');
        hold on;
        Output.DirTotalTime(i) = length(GoodPeakIndices)*TemplateLen/TotalFileTime;
    end
    
    % Next analyse the output data from undirected songs
    if (~isempty(handles(i).UnDirPeaks{SyllIndex}))
        MotifPeaks = handles(i).UnDirPeaks{SyllIndex}(:,1);
        GoodPeakIndices = find(MotifPeaks >= Criterion);
    
        if (~isempty(GoodPeakIndices))
            figure(PlotFigure);
            subplot(3, 1, 1);
            errorbar(i+0.1, mean(MotifPeaks(GoodPeakIndices)), std(MotifPeaks(GoodPeakIndices)), 'bs');
            hold on;
            
            figure(PlotFigure);
            subplot(3, 1, 2);
            GoodPeaks = sort(MotifPeaks(GoodPeakIndices));
            plot(i+0.1, median(GoodPeaks), 'bs');
            hold on;
            plot([i+0.1 i+0.1], [GoodPeaks(ceil(0.25*length(GoodPeaks))) GoodPeaks(ceil(0.75*length(GoodPeaks)))], 'b');
        end
        
        figure(PlotFigure);
        subplot(3, 1, 3);
        TotalFileTime = sum(handles(i).UnDirFileLens);
        plot(i+0.1, length(GoodPeakIndices)*TemplateLen/TotalFileTime, 'bs');
        hold on;
        Output.UnDirTotalTime(i) = length(GoodPeakIndices)*TemplateLen/TotalFileTime;
    end
    
end

figure(PlotFigure);
set(gcf, 'Color', 'w');
subplot(3,1,1);
set(gca, 'Box', 'off');
axis tight;
set(gca, 'FontSize', 12);
ylabel('Template match value', 'FontSize', 12, 'FontWeight', 'bold');

subplot(3,1,2);
set(gca, 'Box', 'off');
axis tight;
set(gca, 'FontSize', 12);
ylabel('Template match value', 'FontSize', 12, 'FontWeight', 'bold');

subplot(3,1,3);
axis tight;
set(gca, 'Box', 'off');
xlabel('Day #', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12);
ylabel('Proportion', 'FontSize', 12, 'FontWeight', 'bold');

disp('Finished analysis');
