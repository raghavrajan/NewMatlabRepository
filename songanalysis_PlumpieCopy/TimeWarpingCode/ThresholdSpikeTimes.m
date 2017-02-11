function [SpikeTimes,SpikeAmplitudes,SpikeWaveforms,ThresholdingParameters] = ThresholdSpikeTimes(DirectoryName,FileNames,ChannelNo,RecordLengths, FileType)

if (strfind(FileType, 'obs'))
    channel_string = ['obs',num2str(ChannelNo),'r'];
end

SpikeTimes = [];
SpikeWaveforms = [];
SpikeAmplitudes = [];

RecTime = 0;

for i = 1:length(FileNames),
    disp([FileNames{i}]);
    if (strfind(FileType, 'obs'))
        [RawData,Fs] = soundin_copy(DirectoryName,FileNames{i},channel_string);
        RawData = RawData * 500/32768; 
    else
        if (strfind(FileType, 'okrank'))
            PresentDirectory = pwd;
            cd(DirectoryName);
            [RawData,Fs] = read_okrank_data(FileNames{i},num2str(ChannelNo));
            RawData = RawData * 100;
            cd(PresentDirectory);
        end
    end
    MeanRawData(i) = mean(RawData);
    STDRawData(i) = std(RawData);
end

if (exist('MeanRawData', 'var'))
    disp(['Mean and std of raw data are ', num2str(max(MeanRawData)), ' and ', num2str(max(STDRawData))]);
end

Prompt = {'Enter the lower threshold (uV)','Enter the upper threshold (uV)','Enter the window size for detecting spikes (no of points)','Enter the window size for skipping detection after a spike (no of points)'};
DialogTitle = 'Input for thresholding parameters';
DefaultParameters = {'-100','100','32','32'};

ThresholdingParameters = inputdlg(Prompt, DialogTitle, 1, DefaultParameters);

LowerThreshold = str2double(ThresholdingParameters{1});
UpperThreshold = str2double(ThresholdingParameters{2});
WindowSize = str2double(ThresholdingParameters{3});
SkipWindowSize = str2double(ThresholdingParameters{4});

for i = 1:length(FileNames),
    disp([FileNames{i}]);
    if (strfind(FileType, 'obs'))
        [RawData,Fs] = soundin_copy(DirectoryName,FileNames{i},channel_string);
        RawData = RawData * 500/32768; 
    else
        if (strfind(FileType, 'okrank'))
            PresentDirectory = pwd;
            cd(DirectoryName);
            [RawData,Fs] = read_okrank_data(FileNames{i},num2str(ChannelNo));
            RawData = RawData * 100;
            cd(PresentDirectory);
        end
    end
    
    StartTime = RecTime;
    % Amplification is 10000 and the acquisition limits are +5 to -5V - so
    % the above line should make the rawdata into units of microVolts.

    disp(['Mean and std of raw data are ', num2str(mean(RawData)), ' and ', num2str(std(RawData))]);
    [TempSpikes, TempSpikeWaveforms] = FindSpikes(RawData,LowerThreshold,UpperThreshold,WindowSize,SkipWindowSize,0,Fs);
    [TempSpikeAmplitudes, TempSpikeWaveforms] = GetSpikeAmplitudes(DirectoryName,FileNames{i},ChannelNo,TempSpikes, FileType);
    SpikeAmplitudes = [SpikeAmplitudes; TempSpikeAmplitudes];
    SpikeWaveforms = [SpikeWaveforms; TempSpikeWaveforms];
    if (size(TempSpikes,1) == 1)
        DataforSaving = [ones(1,length(TempSpikes)); TempSpikes];
    else
        DataforSaving = [ones(length(TempSpikes),1) TempSpikes];
    end
    
    TempSpikes = TempSpikes + RecTime;
    SpikeTimes = [SpikeTimes; TempSpikes];
    RecTime = RecTime + RecordLengths(i);
    
%    save([FileNames{i}(1:(end-5)),'.spk'],'DataforSaving','-ASCII');
end