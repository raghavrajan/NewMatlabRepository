function [] = SeparateMClustSpikeTimes(DirectoryName, FileList, ChanNo)

cd(DirectoryName);

fid = fopen(FileList, 'r');
FileName = fgetl(fid);
RootFileName = FileName;
StartTime = 0;
FileIndex = 0;

while (ischar(FileName(1)))
    FileIndex = FileIndex + 1;
    [RawData1, Fs] = ReadOKrankData(DirectoryName, FileName, 2);
    FileEndTime = length(RawData1)/Fs;
    if (FileIndex == 1)
        SpikeFiles = dir([FileName,'_chan',num2str(ChanNo),'*.times']);
        TotalClusters = length(SpikeFiles);
    end
    
    OutputSpikeTimes = [];
    for i = 1:TotalClusters,
        SpikeTimes = load(SpikeFiles(i).name);
        FileSpikeTimes = SpikeTimes(find((SpikeTimes >= StartTime) & (SpikeTimes <= (StartTime + FileEndTime))));
        FileSpikeTimes = FileSpikeTimes - StartTime;
        if (size(FileSpikeTimes,1) < size(FileSpikeTimes,2))
            FileSpikeTimes = FileSpikeTimes';
        end
        FileSpikeTimes = [(ones(size(FileSpikeTimes))*i) FileSpikeTimes]; 
        OutputSpikeTimes = [OutputSpikeTimes; FileSpikeTimes];
    end
    StartTime = StartTime + (length(RawData1))/Fs;
    
    outputfile = ['Chan',num2str(ChanNo),'_',FileName,'.spk'];
    save(outputfile,'OutputSpikeTimes', '-ASCII');
    FileName = fgetl(fid);
end
