function [] = MakeAllSyllableTemplatesFromFile(DirectoryName, NoteDir, SongFile, FileType, OutputDir, TemplateIncrement, TimeStretch, FreqStretch, Exclusions, varargin)

if (nargin > 9)
    labels = varargin{1};
    onsets = varargin{2};
    offsets = varargin{3};
end

if (nargin > 12)
    FFTWinSize = varargin{4};
    FFTWinOverlap = varargin{5};
    Normalization = varargin{6};
end

% THis is to normalize the whole raw song by rms of the whole song.

if (nargin > 15)
    RawDataNormalization = varargin{7};
end

PresentDir = pwd;
%disp(SongFile);

[Song, Fs] = GetData(DirectoryName, SongFile, FileType, 0);
if (exist('RawDataNormalization', 'var'))
    if (RawDataNormalization == 1)
        Song = Song/sqrt(mean(Song.^2));
    end
end

if (~exist('labels', 'var'))
    cd(DirectoryName);
    cd(NoteDir);
    load([SongFile, '.not.mat']);
end

UniqueLabels = unique(labels);

ExcludedSyllables = zeros(size(UniqueLabels));
for i = 1:length(UniqueLabels),
    if ~isempty(strfind(Exclusions, UniqueLabels(i)))
        ExcludedSyllables(i) = 1;
    end
end

UniqueLabels(find(ExcludedSyllables == 1)) = [];

if (exist('FFTWinSize', 'var'))
    [SyllableTemplates, TemplatePNGFig_Syllables, TemplatePNGFig_SyllLabelTimes] = MakeTemplatesMicrolesionSpectralMatchAnalysis(Song, Fs, [(onsets(:)/1000) (offsets(:)/1000)], labels, UniqueLabels, TimeStretch, FreqStretch, FFTWinSize, FFTWinOverlap, Normalization);
else
    [SyllableTemplates, TemplatePNGFig_Syllables, TemplatePNGFig_SyllLabelTimes] = MakeTemplatesMicrolesionSpectralMatchAnalysis(Song, Fs, [(onsets(:)/1000) (offsets(:)/1000)], labels, UniqueLabels, TimeStretch, FreqStretch);
end

OutputFileName = [SongFile, '.SyllTemplates.', num2str(TemplateIncrement), '.template.mat'];
FileSep = filesep;

save(fullfile(OutputDir, OutputFileName), 'SyllableTemplates');
PlotSpectrogram_SongVar(TemplatePNGFig_Syllables, Fs);
for i = 1:length(UniqueLabels),
    text(TemplatePNGFig_SyllLabelTimes(i), 7000, UniqueLabels(i), 'FontSize', 24, 'Color', 'b');
end
set(gca, 'YColor', 'w');
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(OutputDir, [OutputFileName, '.SyllableTemplates.png']), '-dpng', '-r300');
disp(['Saved templates to ', OutputFileName, ' in directory ', OutputDir]);
cd(PresentDir);