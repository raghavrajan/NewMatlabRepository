function [] = DoAllDaysBoutSpectralMatches(BirdName, RawDataDirs, TemplateDir, TextFileDir, DirTextFiles, UnDirTextFiles, Motif, OutputFileLabel, FileType, StretchVals, PlotOption, DirBoutLimits, UnDirBoutLimits)

PresentDir = pwd;

cd(TemplateDir);

TemplateFiles = dir(['*.', Motif, '.*template.mat']);
clear MotifTemplate;
for j = 1:min([length(TemplateFiles) 3]),
    disp(['Loaded template from ', TemplateFiles(j).name]);
    TempMotifTemplate = load(TemplateFiles(j).name);
    MotifTemplate{j} = TempMotifTemplate.MotifTemplate;
    MotifTemplateRawData{j} = TempMotifTemplate.Motif;
    MotifTemplateRawDataFs{j} = TempMotifTemplate.Fs;
end

if (strfind(FileType, 'wav'))
    GetRidNoise = 'yes';
else
    GetRidNoise = 'no';
end

cd(TextFileDir);
if ((length(DirTextFiles) == length(UnDirTextFiles)) && (length(DirTextFiles) == length(RawDataDirs)))
    for i = 1:length(DirTextFiles),
        DoOneDayBoutSpectralMotifMatches(TextFileDir, MotifTemplate, ['.', Motif, '.', OutputFileLabel], RawDataDirs{i}, DirTextFiles{i}, UnDirTextFiles{i}, FileType, DirBoutLimits{i}, UnDirBoutLimits{i}, StretchVals, GetRidNoise, PlotOption, MotifTemplateRawData, MotifTemplateRawDataFs);
    end
else
    disp('There should be equal number of directed and undirected files and raw data directories');
    return;
end
cd(PresentDir);
disp(['Finished analysis for ', BirdName]);
