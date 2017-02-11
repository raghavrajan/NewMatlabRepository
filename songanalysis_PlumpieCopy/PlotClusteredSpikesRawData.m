function [] = PlotClusteredSpikesRawData(RawFile);

[DirectoryName, FileName, Ext] = fileparts(RawFile);

if (length(DirectoryName) < 1)
    DirectoryName = pwd;
end

if (DirectoryName(end) ~= '/')
    DirectoryName(end + 1) = '/';
end

SpikeTimesFile = dir([FileName,'*.mat']);

[DataCh1, Fs] = soundin_copy(DirectoryName,[FileName, Ext],'obs4r');
DataCh1 = DataCh1 * 1000/32768;

[DataCh2, Fs] = soundin_copy(DirectoryName,[FileName, Ext],'obs3r');
DataCh2 = DataCh2 * 1000/32768;

[DataCh3, Fs] = soundin_copy(DirectoryName,[FileName, Ext],'obs2r');
DataCh3 = DataCh3 * 1000/32768;

[DataCh4, Fs] = soundin_copy(DirectoryName,[FileName, Ext],'obs1r');
DataCh4 = DataCh4 * 1000/32768;

% SpikeTimes = load(SpikeTimesFile.name);

load(SpikeTimesFile.name);
SpikeTimes = T/10000;

Time = 0:1/Fs:((length(DataCh1)-1)/Fs);

SpikeTimeIndices = round(SpikeTimes * Fs);

figure;
subplot(4,1,1);
plot(Time,DataCh1);
subplot(4,1,2);
plot(Time,DataCh2);
subplot(4,1,3);
plot(Time,DataCh3);
subplot(4,1,4);
plot(Time,DataCh4);
hold on;
% for i = 1:length(SpikeTimes),
%     plot(Time((SpikeTimeIndices(i) - 8):(SpikeTimeIndices(i) + 24)),DataCh1((SpikeTimeIndices(i) - 8):(SpikeTimeIndices(i) + 24)),'r-');
% end

