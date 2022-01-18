function [] = ConvertMDAFilesToSpikeTimeFiles(DataDir, FileList, FileType, MDASpikeTimesFile, ChanNoToWrite, ActualChanNo, SpikeClusterNos, TetrodeOrNot)

PresentDir = pwd;
cd(DataDir);

% First read spike times
SpikeTimes = readmda(MDASpikeTimesFile);
% This is a 3 X n array where
% Row 1 is the channel index
% Row 2 is the 

% Now get the file lengths for all the files in the filelist
Fid = fopen(FileList, 'r');
Files = textscan(Fid, '%s', 'DeLimiter', '\n');
Files = Files{1};
fclose(Fid);

for i = 1:length(Files),
    [RawData, Fs] = GetData(DataDir, Files{i}, FileType, 1);
    FileLen(i) = 1000 * length(RawData)/Fs; % in msec
end

CumulativeFileLen = cumsum(FileLen);
CumulativeFileLen = [0; CumulativeFileLen(:)];

% First create spike time directory
DotIndex = find(FileList == '.');
MDAIndex = strfind(FileList, '.mda');
if (isempty(MDAIndex))
    SpikeTimesDir = FileList(1:DotIndex(end)-1);
else
    LastDotBeforeMDA = find(DotIndex < MDAIndex(1), 1, 'last');
    SpikeTimesDir = FileList(1:LastDotBeforeMDA-1);
end

if (~exist(SpikeTimesDir, 'dir'))
    if (exist(SpikeTimesDir, 'file'))
        DotIndex = find(SpikeTimesDir == '.');
        SpikeTimesDir = SpikeTimesDir(1:DotIndex(end)-1);
    end
    mkdir(SpikeTimesDir);
end

% First convert all spiketimes to msec
CumulativeSpikeTimes = SpikeTimes(2,:)/(Fs/1000);

% Now if this is a tetrode, then don't bother about channel no, just pull
% out the spike times for each file and put it into a spike time file
for i = 2:length(CumulativeFileLen),
    if (strcmp(TetrodeOrNot, 'Tetrode'))
        FileSpikeIndices = find((CumulativeSpikeTimes >= CumulativeFileLen(i-1)) & (CumulativeSpikeTimes < CumulativeFileLen(i)));
    else
        FileSpikeIndices = find((CumulativeSpikeTimes >= CumulativeFileLen(i-1)) & (CumulativeSpikeTimes < CumulativeFileLen(i)) & (SpikeTimes(1,:) == ChanNoToWrite));
    end
    FileSpikeTimes = CumulativeSpikeTimes(FileSpikeIndices) - CumulativeFileLen(i-1);
    FileClusterIndices = SpikeTimes(3, FileSpikeIndices);
    
    % Now keep the spike cluster nos that I want to keep and make the rest
    % as 0
    UniqueClusterIndices = unique(FileClusterIndices);
    NoiseClusters = setdiff(UniqueClusterIndices, SpikeClusterNos);
    for j = NoiseClusters(:)',
        FileClusterIndices(find(FileClusterIndices == j)) = 0;
    end
    
    % Now write it to a spike time file
    Fid = fopen(fullfile(SpikeTimesDir, ['Chan', num2str(ActualChanNo), '_', Files{i-1}, '.spk']), 'w');
    for j = 1:length(FileSpikeTimes),
        fprintf(Fid, '%g\t%g\n', FileClusterIndices(j), FileSpikeTimes(j));
    end
    fclose(Fid);
end
cd(PresentDir);