function [CorrectMatches] = MA_PlotMatchSpectrograms(Parameters, Indices, i, SyllTemp, TemplateMatchValues, TemplateLen, PreOrPost, Context)

NumRows = 5;
NumCols = 5;

PlotIndices = [];

figure;
set(gcf, 'Position', [108 111 1212 578]);
subplot(NumRows, NumCols, NumCols + 1);
PlotIndices = [PlotIndices (NumCols + 1)];

SyllableData = Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).RawSound;
SyllableDataFs = Parameters.SyllableTemplate.SyllableTemplates{SyllTemp}{1}.MotifTemplate(1).Fs;
PlotSpectrogramInAxis_SongVar(SyllableData, (1:1:length(SyllableData))/SyllableDataFs, SyllableDataFs, gca);
Temp = axis;
PlotDur(1,:) = Temp(2) - Temp(1);

title('Syllable template', 'FontSize', 8, 'FontWeight', 'bold');

RandSubset = randperm(length(Indices));

Index = 1;
for k = 1:min([20 length(RandSubset)]),
    if (Index > 5)
        subplot(NumRows, NumCols, Index+NumCols);
        PlotIndices = [PlotIndices (Index + NumCols)];
    else
        subplot(NumRows, NumCols,Index);
        PlotIndices = [PlotIndices Index];
    end
    FileIndex = TemplateMatchValues(Indices(RandSubset(k)), 3);
    
    DataDir = [eval(['Parameters.', PreOrPost, 'DataDir{i}']), '/'];
    SongFileName = eval(['Parameters.', PreOrPost, Context, 'SongFileNames{i}{FileIndex}']);
    
    disp([SongFileName, ': Match value ', num2str(TemplateMatchValues(Indices(RandSubset(k)),1)), '@ ', num2str(TemplateMatchValues(Indices(RandSubset(k)), 2))]);
    
    Onsets = eval(['Parameters.', PreOrPost, Context, 'Onsets{i}{FileIndex}']);
    Offsets = eval(['Parameters.', PreOrPost, Context, 'Offsets{i}{FileIndex}']);
    BoutLen = eval(['Parameters.', PreOrPost, Context, 'Lens{i}{FileIndex}']);
    [OnsetTime, OffsetTime] = MA_FindMatchOnsetsOffsets(Onsets, Offsets, TemplateMatchValues(Indices(RandSubset(k)),:), TemplateLen, BoutLen);
    
    PlotSpectrogramInAxis(DataDir, SongFileName, Parameters.FileType, gca, [OnsetTime OffsetTime]);
    hold on;
    plot(ones(1,2) * (TemplateMatchValues(Indices(RandSubset(k)),2)), [300 8000], 'b--');
    title([SongFileName, ': ', num2str(TemplateMatchValues(Indices(RandSubset(k)),2))], 'FontSize', 8);
    Temp = axis;
    PlotDur(end+1) = Temp(2) - Temp(1);
    Index = Index + 1;
end

for k = PlotIndices(:)',
    subplot(NumRows, NumCols, k);
    Temp = axis;
    Temp(2) = Temp(1) + max(PlotDur);
    axis(Temp);
end
    
CorrectMatches = inputdlg('Enter the number of plots that have correct matches', 'Correct plot selection box');
