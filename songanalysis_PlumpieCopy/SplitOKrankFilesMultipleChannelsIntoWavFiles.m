function [] = SplitOKrankFilesMultipleChannelsIntoWavFiles(DirectoryName, FileList)

Fid = fopen(FileList, 'r');
FileNames = textscan(Fid, '%s', 'DeLimiter', '\n');
FileNames = FileNames{1};
fclose(Fid);

for i = 1:length(FileNames),
    for j = 1:2,
        [RawData, Fs] = GetData(DirectoryName, FileNames{i}, 'okrank', j);
    
        NewFs = 44100;
        RawData = RawData/10;
        RawData = resample(RawData,44100,Fs);
    
        OutputFileName = [FileNames{i}, '.Chan', num2str(j), '.wav'];
        audiowrite(OutputFileName, RawData, NewFs);
    end
end

disp(['Split ', num2str(length(FileNames)), ' okrank files']);
        
        
    
    
    