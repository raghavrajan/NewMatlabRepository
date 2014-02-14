function [] = PrepareObsFilesForMClust(DirectoryName,RawDataFile,UpperThreshold,WindowSize,WindowToSkip)

PresentDirectory = pwd;

if (DirectoryName(end) ~= '/')
    DirectoryName(end + 1) = '/';
end

cd(DirectoryName);

[DataCh1, Fs] = soundin_copy(DirectoryName,RawDataFile,'obs4r');
DataCh1 = DataCh1 * 1000/32768;

[DataCh2, Fs] = soundin_copy(DirectoryName,RawDataFile,'obs3r');
DataCh2 = DataCh2 * 1000/32768;

[DataCh3, Fs] = soundin_copy(DirectoryName,RawDataFile,'obs2r');
DataCh3 = DataCh3 * 1000/32768;

[DataCh4, Fs] = soundin_copy(DirectoryName,RawDataFile,'obs1r');
DataCh4 = DataCh4 * 1000/32768;

StartTime = 0;
[ST, wv1, wv2, wv3, wv4] = tetrode_find_spikes(DataCh1,DataCh2,DataCh3,DataCh4,UpperThreshold,WindowSize,StartTime,WindowToSkip);

wv1 = wv1';
wv2 = wv2';
wv3 = wv3';
wv4 = wv4';

% Times = load('TempSpikeTimes.txt');
% Ch1Waveforms = load('TempSpikeWaveformsCh1.txt');
% Ch2Waveforms = load('TempSpikeWaveformsCh2.txt');
% Ch3Waveforms = load('TempSpikeWaveformsCh3.txt');
% Ch4Waveforms = load('TempSpikeWaveformsCh4.txt');

TotalNoofSpikes = length(ST);

T = ST;
T = T * 10000;
wv = [];

SpikeCount = 1:1:TotalNoofSpikes;
SpikeCount = SpikeCount';

wv(SpikeCount,1,:) = wv1;
wv(SpikeCount,2,:) = wv2;
wv(SpikeCount,3,:) = wv3;
wv(SpikeCount,4,:) = wv4;

outputfile = strcat(DirectoryName,RawDataFile,'.mat');

save(outputfile,'T','wv');

disp('Finished');