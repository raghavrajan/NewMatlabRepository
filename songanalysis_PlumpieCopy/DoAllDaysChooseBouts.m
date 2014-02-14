function [DirBoutLimits, UnDirBoutLimits] = DoAllDaysChooseBouts(BirdName, RawDataDirs, TextFileDir, DirTextFiles, UnDirTextFiles, FileType, ChooseBouts)

PresentDir = pwd;

if (strfind(FileType, 'wav'))
    GetRidNoise = 'yes';
else
    GetRidNoise = 'no';
end

cd(TextFileDir);
if ((length(DirTextFiles) == length(UnDirTextFiles)) && (length(DirTextFiles) == length(RawDataDirs)))
    for i = 1:length(DirTextFiles),
        cd(TextFileDir);
        [DirBoutLimits{i}] = BoutMotifSpectralMatch_ChooseBouts(RawDataDirs{i}, DirTextFiles{i}, FileType, ChooseBouts, GetRidNoise);
    end
    for i = 1:length(UnDirTextFiles),
        cd(TextFileDir);
        [UnDirBoutLimits{i}] = BoutMotifSpectralMatch_ChooseBouts(RawDataDirs{i}, UnDirTextFiles{i}, FileType, ChooseBouts, GetRidNoise);
    end
else
    disp('There should be equal number of directed and undirected files and raw data directories');
    return;
end
cd(PresentDir);
disp(['Finished choosing bouts for ', BirdName]);
