function [] = PlotContinuousNeuronFiringRateWithContext(FileList, FileType, SpikeFileDir, ChanNo, ClusterNos, DirectedTimes, PlotFig, StartTime, DT)

PresentDir = pwd;

% First read in all the filenames

Fid = fopen(FileList, 'r');
TempFiles = textscan(Fid, '%s', 'DeLimiter', '\n');
TempFiles = TempFiles{1};
fclose(Fid);

% Now load up spike times and get total file length

SpikeTimes = [];
for i = 1:length(TempFiles),
    cd(SpikeFileDir);
    SpikeInfo = load(['Chan', num2str(ChanNo), '_', TempFiles{i}, '.spk']);
    TempSpikeTimes = [];
    for j = 1:length(ClusterNos),
        TempSpikeTimes = [TempSpikeTimes; SpikeInfo(find(SpikeInfo(:,1) == ClusterNos(j)), 2)];
    end
    SpikeTimes = [SpikeTimes; (TempSpikeTimes(:)*1000 + StartTime)];
    
    cd(PresentDir);
    
    [Data, Fs] = ASSLGetRawData(PresentDir, TempFiles{i}, FileType, ChanNo);
    DataLen = length(Data)*1000/Fs ;
    StartTime = StartTime + DataLen;
end


% Now calculate instantaneous firing rate
IFR = 1./diff(SpikeTimes);
IFR = IFR(:);
figure(PlotFig);
plot(SpikeTimes(2:end), IFR, 'k');