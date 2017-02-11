function [] = PlotOtherSyllableMatchValues(handles, Syllable, varargin)


% File to plot the spectrograms of template match values for each of the
% days

BinSize = 0.01;

if (nargin > 2)
    InputData = varargin{1};
end

if (nargin > 3)
    Syllables = varargin{2};
end

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

Temp = inputdlg('Enter the threshold above which peaks should be considered (0.5 is a good place to start)', 'Peak detection threshold');

MPH = str2double(Temp{1});

Colours = 'rgbcmky';

SyllIndex = intersect(find(cellfun(@length, strfind([handles(1).Labels], Syllable))), find(cellfun(@length, handles(1).Labels) == length(Syllable)));

for i = 1:length(Syllables),
    if (i ~= SyllIndex)
        PlotFigure(i) = figure;
        set(PlotFigure(i), 'Position', [427 170 675 525]);
    end
end

for i = 1:length(handles),
    Legend{i} = ['Day #', num2str(i)];
 
    % First analyse the output data from directed songs
    
    if (~isempty(handles(i).DirPeaks{SyllIndex}))
        MotifPeakIndices = find(handles(i).DirPeaks{SyllIndex}(:,1) >= MPH);
        for k = 1:length(Syllables),
            PeaksAndTimes = [];
            PrePeaks = [];
            PostPeaks = [];
            
            if (k ~= SyllIndex)
                for j = 1:length(MotifPeakIndices),
                    MatchTime = handles(i).DirPeaks{SyllIndex}(MotifPeakIndices(j),2);
        
                    SyllIndices = find(handles(i).DirPeaks{k}(:,3) == handles(i).DirPeaks{SyllIndex}(MotifPeakIndices(j),3));
                    MatchIndices = SyllIndices(find((handles(i).DirPeaks{k}(SyllIndices,2) >= (MatchTime - 0.5)) & (handles(i).DirPeaks{k}(SyllIndices,2) <= (MatchTime + TemplateLen + 0.5))));
                    PeaksAndTimes = [PeaksAndTimes; [(handles(i).DirPeaks{k}(MatchIndices,2) - MatchTime) (handles(i).DirPeaks{k}(MatchIndices,1))]];
                    
                    MatchIndices = SyllIndices(find((handles(i).DirPeaks{k}(SyllIndices,2) >= (MatchTime - 0.5)) & (handles(i).DirPeaks{k}(SyllIndices,2) <= (MatchTime))));
                    if (~isempty(MatchIndices))
                        PrePeaks = [PrePeaks; [handles(i).DirPeaks{SyllIndex}(MotifPeakIndices(j),1) max(handles(i).DirPeaks{k}(MatchIndices,1))]];
                    end                        
                    
                    MatchIndices = SyllIndices(find((handles(i).DirPeaks{k}(SyllIndices,2) >= (MatchTime)) & (handles(i).DirPeaks{k}(SyllIndices,2) <= (MatchTime + TemplateLen + 0.5))));
                    if (~isempty(MatchIndices))
                        PostPeaks = [PostPeaks; [handles(i).DirPeaks{SyllIndex}(MotifPeakIndices(j),1) max(handles(i).DirPeaks{k}(MatchIndices,1))]];
                    end
                end
            
                figure(PlotFigure(k));
                
                subplot(2,3,1);
                hold on;
                
                Edges = 0:BinSize:max(PeaksAndTimes);
                MeanPeaksAndTimes = [];
                for Edge = 1:length(Edges)-1,
                    MeanPeaksAndTimes(Edge,1) = Edges(Edge) + BinSize/2;
                    MeanPeaksAndTimes(Edge,2) = mean(PeaksAndTimes(find((PeaksAndTimes(:,1) >= Edges(Edge)) & (PeaksAndTimes(:,1) < Edges(Edge+1))),2));
                    MeanPeaksAndTimes(Edge,3) = std(PeaksAndTimes(find((PeaksAndTimes(:,1) >= Edges(Edge)) & (PeaksAndTimes(:,1) < Edges(Edge+1))),2));
                end
                plot(MeanPeaksAndTimes(:,1), MeanPeaksAndTimes(:,2), [Colours(i), '-o']);
                
                subplot(2,3,2);
                hold on;
                plot(PrePeaks(:,1), PrePeaks(:,2), [Colours(i), 'o']);
                
                subplot(2,3,3);
                hold on;
                plot(PostPeaks(:,1), PostPeaks(:,2), [Colours(i), 'o']);
            end
        end
    end
    
    % Next analyse the output data from undirected songs
    if (~isempty(handles(i).UnDirPeaks{SyllIndex}))
        MotifPeakIndices = find(handles(i).UnDirPeaks{SyllIndex}(:,1) >= MPH);
        for k = 1:length(Syllables),
            PeaksAndTimes = [];
            PrePeaks = [];
            PostPeaks = [];
            
            if (k ~= SyllIndex)
                for j = 1:length(MotifPeakIndices),
                    MatchTime = handles(i).UnDirPeaks{SyllIndex}(MotifPeakIndices(j),2);
        
                    SyllIndices = find(handles(i).UnDirPeaks{k}(:,3) == handles(i).UnDirPeaks{SyllIndex}(MotifPeakIndices(j),3));
                    MatchIndices = SyllIndices(find((handles(i).UnDirPeaks{k}(SyllIndices,2) >= (MatchTime - 0.5)) & (handles(i).UnDirPeaks{k}(SyllIndices,2) <= (MatchTime + TemplateLen + 0.5))));
                    PeaksAndTimes = [PeaksAndTimes; [(handles(i).UnDirPeaks{k}(MatchIndices,2) - MatchTime) (handles(i).UnDirPeaks{k}(MatchIndices,1))]];
                    
                    MatchIndices = SyllIndices(find((handles(i).UnDirPeaks{k}(SyllIndices,2) >= (MatchTime - 0.5)) & (handles(i).UnDirPeaks{k}(SyllIndices,2) <= (MatchTime))));
                    if (~isempty(MatchIndices))
                        PrePeaks = [PrePeaks; [handles(i).UnDirPeaks{SyllIndex}(MotifPeakIndices(j),1) max(handles(i).UnDirPeaks{k}(MatchIndices,1))]];
                    end                
                    MatchIndices = SyllIndices(find((handles(i).UnDirPeaks{k}(SyllIndices,2) >= (MatchTime)) & (handles(i).UnDirPeaks{k}(SyllIndices,2) <= (MatchTime + TemplateLen + 0.5))));
                    if (~isempty(MatchIndices))
                        PostPeaks = [PostPeaks; [handles(i).UnDirPeaks{SyllIndex}(MotifPeakIndices(j),1) max(handles(i).UnDirPeaks{k}(MatchIndices,1))]];
                    end
                    
                end
            
                figure(PlotFigure(k));
                
                subplot(2,3,4);
                Edges = 0:BinSize:max(PeaksAndTimes);
                MeanPeaksAndTimes = [];
                for Edge = 1:length(Edges)-1,
                    MeanPeaksAndTimes(Edge,1) = Edges(Edge) + BinSize/2;
                    MeanPeaksAndTimes(Edge,2) = mean(PeaksAndTimes(find((PeaksAndTimes(:,1) >= Edges(Edge)) & (PeaksAndTimes(:,1) < Edges(Edge+1))),2));
                    MeanPeaksAndTimes(Edge,3) = std(PeaksAndTimes(find((PeaksAndTimes(:,1) >= Edges(Edge)) & (PeaksAndTimes(:,1) < Edges(Edge+1))),2));
                end
                hold on;
                plot(MeanPeaksAndTimes(:,1), MeanPeaksAndTimes(:,2), [Colours(i), '-o']);
                                
                subplot(2,3,5);
                hold on;
                plot(PrePeaks(:,1), PrePeaks(:,2), [Colours(i), 'o']);
                
                subplot(2,3,6);
                hold on;
                plot(PostPeaks(:,1), PostPeaks(:,2), [Colours(i), 'o']);
                
            end
        end
    end
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
    
Flag = 1;

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