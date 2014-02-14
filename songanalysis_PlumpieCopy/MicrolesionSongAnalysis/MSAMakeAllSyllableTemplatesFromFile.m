function [Template, TemplateLabels] = MSAMakeAllSyllableTemplatesFromFile(DirectoryName, NoteDir, SongFile, FileType, OutputDir, TemplateIncrement)

PresentDir = pwd;
cd(DirectoryName);
disp(SongFile);

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

cd(NoteDir);
load([SongFile, '.not.mat']);
UniqueLabels = unique(labels);

for i = 1:length(UniqueLabels),
    Indices = find(labels == UniqueLabels(i));
    for j = 1:length(Indices),
        Template{i}{j} = MakeTemplatesMicrolesionSpectralMatchAnalysis(Song, Fs, [(onsets(Indices(j))/1000 - 0.01) (offsets(Indices(j))/1000 + 0.01)], OutputDir, [SongFile, '.', labels(Indices(j)), '.', num2str(j + TemplateIncrement), '.template.mat']);
        TemplateLabels(i) = UniqueLabels(i);
        title([labels(Indices(j)), num2str(j)], 'FontSize', 16, 'FontWeight', 'bold');
    end
end