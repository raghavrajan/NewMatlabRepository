function [] = PlotMatchValDistributions(handles, Syllable, varargin)

% File to plot the distribution of template match values for each of the
% days

Temp = inputdlg('Enter the threshold above which peaks should be considered (0.5 is a good place to start)', 'Peak detection threshold');

MPH = str2double(Temp{1});

Colours = 'rgbcmky';

PlotFigure = figure;


SyllIndex = intersect(find(cellfun(@length, strfind([handles(1).Labels], Syllable))), find(cellfun(@length, handles(1).Labels) == length(Syllable)));
for i = 1:length(handles),
    Legend{i} = ['Day #', num2str(i)];
 
    % First analyse the output data from directed songs

    if (isfield(handles(i), 'DirPeaks'))
        MotifPeaks = handles(i).DirPeaks{SyllIndex}(:,1);
        MotifPeaks = MotifPeaks(find(MotifPeaks >= MPH));
        figure(PlotFigure);
        subplot(2, 1, 1);
        plot(0:0.1:max(MotifPeaks), 100*histc(MotifPeaks, (0:0.1:max(MotifPeaks)))/sum(histc(MotifPeaks, (0:0.1:max(MotifPeaks)))), Colours(i), 'LineWidth', 2);
        hold on;
    end
    
    % Next analyse the output data from undirected songs
    MotifPeaks = handles(i).UnDirPeaks{SyllIndex}(:,1);
    MotifPeaks = MotifPeaks(find(MotifPeaks >= MPH));
    figure(PlotFigure);
    subplot(2, 1, 2);
    plot(0:0.1:max(MotifPeaks), 100*histc(MotifPeaks, (0:0.1:max(MotifPeaks)))/sum(histc(MotifPeaks, (0:0.1:max(MotifPeaks)))), Colours(i), 'LineWidth', 2);
    hold on;
end

% Now for all the random comparisons
MotifPeaks = handles(1).RandomSongComparisonUnDirPeaks{SyllIndex}(:,1);
MotifPeaks = MotifPeaks(find(MotifPeaks >= MPH));
figure(PlotFigure);
subplot(2, 1, 2);
plot(0:0.1:max(MotifPeaks), 100*histc(MotifPeaks, (0:0.1:max(MotifPeaks)))/sum(histc(MotifPeaks, (0:0.1:max(MotifPeaks)))), [Colours(1), ':'], 'LineWidth', 2);
hold on;

MotifPeaks = handles(1).RandomSoundComparisonUnDirPeaks{SyllIndex}(:,1);
MotifPeaks = MotifPeaks(find(MotifPeaks >= MPH));
figure(PlotFigure);
subplot(2, 1, 2);
plot(0:0.1:max(MotifPeaks), 100*histc(MotifPeaks, (0:0.1:max(MotifPeaks)))/sum(histc(MotifPeaks, (0:0.1:max(MotifPeaks)))), [Colours(1), '--'], 'LineWidth', 2);
hold on;

MotifPeaks = handles(1).ShuffledSongComparisonsUnDirPeaks{SyllIndex}(:,1);
MotifPeaks = MotifPeaks(find(MotifPeaks >= MPH));
figure(PlotFigure);
subplot(2, 1, 2);
plot(0:0.1:max(MotifPeaks), 100*histc(MotifPeaks, (0:0.1:max(MotifPeaks)))/sum(histc(MotifPeaks, (0:0.1:max(MotifPeaks)))), [Colours(1), '-.'], 'LineWidth', 2);
hold on;

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
axis([0 max(TempAxis(:,2))*1.05 0 max(TempAxis(:,4))*1.05]);
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
axis([0 max(TempAxis(:,2))*1.05 0 max(TempAxis(:,4))*1.05]);
if (length(Syllable) > 1)
    title(['UnDirected motif ', Syllable], 'FontSize', 12, 'FontWeight', 'bold');
else
    title(['UnDirected syllable ', Syllable], 'FontSize', 12, 'FontWeight', 'bold');
end
xlabel('Template match value', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('%', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12);
    
disp('Finished analysis');