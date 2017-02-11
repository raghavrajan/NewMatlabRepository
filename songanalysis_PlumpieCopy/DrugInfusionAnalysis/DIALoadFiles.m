function [DataFiles] = DIALoadFiles(DIAData)

BirdName = DIAData.BirdName;
ExpDate = DIAData.ExpDate;
ExpStartTime = str2double(DIAData.ExpStartTime);
ExpEndTime = str2double(DIAData.ExpEndTime);

Files = dir([BirdName, '_', ExpDate, '*.rec']);

DataFiles = [];

for i = 1:length(Files),
    DataFileName = Files(i).name(1:(end -4));
    DataFileTime = str2double(DataFileName((end - 5):end));
    if ((DataFileTime > ExpStartTime) && (DataFileTime < ExpEndTime))
        DataFiles = [DataFiles; DataFileName];
    end
end

disp(['Loaded ', num2str(size(DataFiles, 1)), ' files']);