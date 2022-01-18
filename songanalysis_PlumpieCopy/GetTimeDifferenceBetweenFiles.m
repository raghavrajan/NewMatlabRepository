function [FileTimeDifferences] = GetTimeDifferenceBetweenFiles(DataDir, BaseFileName, FileType, FemaleIntroductionTime)

%=========================================================================
% This is a function to get the time difference between all the files in a
% directory. Time of a file is taken as the last 6 digits of the name
% before the extension
%=========================================================================

if (length(FemaleIntroductionTime) == 5)
    FemaleIntroTime = 3600*str2double(FemaleIntroductionTime(end-4)) + 60*str2double(FemaleIntroductionTime(end-3:end-2)) + str2double(FemaleIntroductionTime(end-1:end));
else
    FemaleIntroTime = 3600*str2double(FemaleIntroductionTime(end-5:end-4)) + 60*str2double(FemaleIntroductionTime(end-3:end-2)) + str2double(FemaleIntroductionTime(end-1:end));
end

FileTime = [];
PresentDir = pwd;
cd(DataDir);
Files = dir([BaseFileName(1:3), '*.', FileType]);
cd(PresentDir);

for i = 1:length(Files),
    switch FileType
        case 'wav'
            FileTime(i) = 3600*str2double(Files(i).name(end-9:end-8)) + 60*str2double(Files(i).name(end-7:end-6)) + str2double(Files(i).name(end-5:end-4));
    end
end

if (isempty(FileTime))
    FileTimeDifferences = [];
else
    SortedFileTimes = sort(FileTime);
    SortedFileTimes(end) = [];
    SortedFileTimes = SortedFileTimes(find(SortedFileTimes > FemaleIntroTime));
    FileTimeDifferences = diff(SortedFileTimes) - 30;
end