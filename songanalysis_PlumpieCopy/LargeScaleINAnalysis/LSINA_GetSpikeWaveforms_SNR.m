function [SpikeWaveforms, SNR] = LSINA_GetSpikeWaveforms_SNR(SURecordingDetails)

PreNSamples = 20;
PostNSamples = 36;

SpikeWaveforms = [];
AlignedSpikeWaveforms = [];
SNR = [];

% First load up all the waveforms
for i = 1:length(SURecordingDetails.ClusterSpikeTimes),
    if (~isempty(SURecordingDetails.ClusterSpikeTimes{i}))
        if strfind(SURecordingDetails.FileType, 'obs')
            [RawSpikeData, Fs] = GetData(SURecordingDetails.DataDirectory, SURecordingDetails.SongFileNames{i}, SURecordingDetails.FileType, SURecordingDetails.SpikeChanNo);
        else
            [RawSpikeData, Fs] = GetData(SURecordingDetails.DataDirectory, SURecordingDetails.SongFileNames{i}, SURecordingDetails.FileType, SURecordingDetails.SpikeChanNo);
        end
        
        % Filter the rawdata between 600 and 6000 Hz - as in some birds I
        % have recorded with the filters open for LFP too.
        
        [b, a] = butter(4, [600*2/Fs 6000*2/Fs]);
        RawSpikeData = filtfilt(b, a, RawSpikeData);
        
        RawSpikeData = RawSpikeData(:)';
        % Now use the inverted or not to invert the data
        if (SURecordingDetails.Spikeinvertedornot == 1)
            RawSpikeData = -RawSpikeData;
        end
        
        % Now adjust for the gain
        RawSpikeData = RawSpikeData/SURecordingDetails.Gain;
        
        TempSpikeTimeIndices = round(SURecordingDetails.ClusterSpikeTimes{i}(:) * Fs/1000);
        % For now I'm just going to get rid of spikes that don't have
        % enough pre and post data - this should happen only for continuous
        % files where typically I would get this data from the previous
        % file or the next file.
        TempSpikeTimeIndices(find(TempSpikeTimeIndices <= PreNSamples)) = [];
        TempSpikeTimeIndices(find(TempSpikeTimeIndices > (length(RawSpikeData) - PostNSamples))) = [];
        
        clear SpikeTimeIndices;
        for j = -PreNSamples:1:PostNSamples,
            SpikeTimeIndices(:,j + PreNSamples + 1) = TempSpikeTimeIndices(:) + j;
        end
        SpikeWaveforms = [SpikeWaveforms; RawSpikeData(SpikeTimeIndices)];
    end
end

% Now to align the waveforms
MaxShift = 5;
Threshold = 0;
FinalPreNSamples = 16;
FinalPostNSamples = 32;

[AlignedSpikeWaveforms] = AlignSpikeWaveformsByMax_SamplesInputArguments(SpikeWaveforms, MaxShift, Threshold, PreNSamples, FinalPreNSamples, FinalPostNSamples);
MeanWF = mean(AlignedSpikeWaveforms);
MeanSubtractedWF = AlignedSpikeWaveforms - repmat(MeanWF, size(AlignedSpikeWaveforms, 1), 1);
SNR = [((max(MeanWF) - min(MeanWF))/std(MeanSubtractedWF(:)*2)) ((max(MeanWF))/std(MeanSubtractedWF(:))) ((min(MeanWF))/std(MeanSubtractedWF(:)))];
