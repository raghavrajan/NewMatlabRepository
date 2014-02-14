function [] = MakeAllSyllableTemplatesFromFile(DirectoryName, NoteDir, SongFile, FileType, OutputDir, TemplateIncrement, TimeStretch, FreqStretch, Exclusions, varargin)

if (nargin > 9)
    labels = varargin{1};
    onsets = varargin{2};
    offsets = varargin{3};
end

PresentDir = pwd;
cd(DirectoryName);
%disp(SongFile);

if (strfind(FileType, 'okrank'))
    [Song, Fs] = ReadOKrankData(DirectoryName, SongFile, 1);
else
    if (strfind(FileType, 'obs'))
        if ispc
            [Song, Fs] = soundin_copy([DirectoryName, '\'], SongFile, 'obs0r');
        else
            [Song, Fs] = soundin_copy([DirectoryName, '/'], SongFile, 'obs0r');
        end
        Song = Song*5/32768;
    else
        if (strfind(FileType, 'wav'))
            cd(DirectoryName);
            [Song, Fs] = wavread(SongFile);
            cd(PresentDir);
        end
    end
end

if (~exist('labels', 'var'))
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

for i = 1:length(UniqueLabels),
    Indices = find(labels == UniqueLabels(i));
    for j = 1:length(Indices),
        SyllableTemplates{i}{j}.MotifTemplate = MakeTemplatesMicrolesionSpectralMatchAnalysis(Song, Fs, [(onsets(Indices(j))/1000) (offsets(Indices(j))/1000)], UniqueLabels(i), TimeStretch, FreqStretch);
        %title([labels(Indices(j)), num2str(j)], 'FontSize', 16, 'FontWeight', 'bold');
    end
end

OutputFileName = [SongFile, '.SyllTemplates.', num2str(TemplateIncrement), '.template.mat'];
FileSep = filesep;

save([OutputDir, FileSep, OutputFileName], 'SyllableTemplates');
disp(['Saved templates to ', OutputFileName, ' in directory ', OutputDir]);