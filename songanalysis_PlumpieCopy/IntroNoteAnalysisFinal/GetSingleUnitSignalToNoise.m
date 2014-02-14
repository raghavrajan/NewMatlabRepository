function [SignalToNoise, AlignedWaveforms, MeanSpikeWaveform, Fs] = GetSingleUnitSignalToNoise(Neural_INR, DataDirectory, SpikeDir, SpikeChanNo, ClusterNos, Threshold, FileType, varargin)


if (nargin > 7)
    if (~isempty(varargin{1}))
        BoutGain = varargin{1};
    end
    if (nargin > 8)
        if (~isempty(varargin{2}))
            BoutSpikeChanNo = varargin{2};
        end
    end
    if (nargin > 9)
        if (~isempty(varargin{3}))
            InvertWaveform = 1;
        end
    end
end

if (~exist('InvertWaveform', 'var'))
    InvertWaveform = 0;
end
   
c;

SpikeIndex = 0;
PreSamples = 24;
PostSamples = 48;

WFIndices = -PreSamples:1:(PostSamples-1);
WFIndices = WFIndices(:)';

Waveforms = [];

cd(DataDirectory);
FileIndex = 0;
FileList = [];
if (exist('BoutGain', 'var'))
	FileBoutGain = [];
end
if (exist('BoutSpikeChanNo', 'var'))
	FileSpikeChanNo = [];
end
for i = 1:length(Neural_INR.BoutDetails)
    if (iscell(Neural_INR.BoutDetails(i).SongFile))
        if (~isempty(FileList))
            if (isempty(find(cellfun(@length, strfind(FileList, Neural_INR.BoutDetails(i).SongFile{1})))))
                FileList = [FileList; Neural_INR.BoutDetails(i).SongFile];
                if (exist('BoutGain', 'var'))
                    FileBoutGain = [FileBoutGain; ones(length(Neural_INR.BoutDetails(i).SongFile), 1)*BoutGain(i)];
                end
                if (exist('BoutSpikeChanNo', 'var'))
                    FileSpikeChanNo = [FileSpikeChanNo; ones(length(Neural_INR.BoutDetails(i).SongFile), 1)*BoutSpikeChanNo(i)];
                end
                SpikeRawData = [];
            end
        else
            SpikeRawData = [];
            FileList = [FileList; Neural_INR.BoutDetails(i).SongFile];
            if (exist('BoutGain', 'var'))
                FileBoutGain = [FileBoutGain; ones(length(Neural_INR.BoutDetails(i).SongFile), 1)*BoutGain(i)];
            end
            if (exist('BoutSpikeChanNo', 'var'))
                FileSpikeChanNo = [FileSpikeChanNo; ones(length(Neural_INR.BoutDetails(i).SongFile), 1)*BoutSpikeChanNo(i)];
            end
        end
        
        if (isempty(SpikeRawData))
            for j = 1:length(Neural_INR.BoutDetails(i).SongFile),
                if (strfind(FileType, 'okrank'))
                    if (exist('BoutSpikeChanNo', 'var'))
                        [TempData, Fs] = SSAReadOKrankData(DataDirectory, DataDirectory, Neural_INR.BoutDetails(i).SongFile{j}, BoutSpikeChanNo(i));
                    else
                        [TempData, Fs] = SSAReadOKrankData(DataDirectory, DataDirectory, Neural_INR.BoutDetails(i).SongFile{j}, SpikeChanNo);
                    end
                else
                   if (strfind(FileType, 'obs'))
                        channel_string = strcat('obs',num2str(SpikeChanNo),'r');
                        [TempData, Fs] = SSASoundIn([DataDirectory, '/'], DataDirectory, Neural_INR.BoutDetails(i).SongFile{j}, channel_string);
                        % Convert to V - 5V on the data acquisition is 32768
                        TempData = TempData * 5/32768;
                    end
                end
                [b, a] = butter(4, [600/16000 6000/16000]);
                TempData = filtfilt(b, a, TempData);         
                TempData = TempData(:);
                if (exist('BoutGain', 'var'))
                    TempData = TempData*BoutGain(i);
                end
                SpikeRawData = [SpikeRawData; TempData];
            end
        end
        
        SpikeIndices = find((Neural_INR.BoutDetails(i).SpikeTimes >= (Neural_INR.BoutDetails(i).onsets(1) - 0.05)) & (Neural_INR.BoutDetails(i).SpikeTimes >= (Neural_INR.BoutDetails(i).offsets(1))));
        SpikeIndices = round(Neural_INR.BoutDetails(i).SpikeTimes(SpikeIndices) * Fs);
        SpikeIndices = SpikeIndices(:);
        SpikeWaveformIndices = repmat(SpikeIndices, 1, length(WFIndices));
        SpikeWaveformIndices = SpikeWaveformIndices + repmat(WFIndices, size(SpikeWaveformIndices,1), 1);
        if (size(SpikeWaveformIndices,1) == 1)
            Waveforms = [Waveforms; SpikeRawData(SpikeWaveformIndices)'];
        else
            Waveforms = [Waveforms; SpikeRawData(SpikeWaveformIndices)];
        end
    else
        SpikeRawData = [];
        if (~isempty(FileList))
            if (iscell(FileList))
                if (isempty(strfind(FileList{end}, Neural_INR.BoutDetails(i).SongFile)))
                    FileList = [FileList; Neural_INR.BoutDetails(i).SongFile];
                    if (exist('BoutGain', 'var'))
                        FileBoutGain = [FileBoutGain; BoutGain(i)];
                    end
                end
                if (exist('BoutSpikeChanNo', 'var'))
                    FileSpikeChanNo = [FileSpikeChanNo; ones(length(Neural_INR.BoutDetails(i).SongFile), 1)*BoutSpikeChanNo(i)];
                end
            else
                if (isempty(strfind(FileList(end,:), Neural_INR.BoutDetails(i).SongFile)))
                    FileList = [FileList; Neural_INR.BoutDetails(i).SongFile];
                    if (exist('BoutGain', 'var'))
                        FileBoutGain = [FileBoutGain; BoutGain(i)];
                    end
                end
                if (exist('BoutSpikeChanNo', 'var'))
                    FileSpikeChanNo = [FileSpikeChanNo; ones(length(Neural_INR.BoutDetails(i).SongFile), 1)*BoutSpikeChanNo(i)];
                end
            end
        else
            FileList = [FileList; Neural_INR.BoutDetails(i).SongFile];
            if (exist('BoutGain', 'var'))
                FileBoutGain = [FileBoutGain; BoutGain(i)];
            end
            if (exist('BoutSpikeChanNo', 'var'))
                FileSpikeChanNo = [FileSpikeChanNo; ones(length(Neural_INR.BoutDetails(i).SongFile), 1)*BoutSpikeChanNo(i)];
            end
        end
        if (strfind(FileType, 'okrank'))
            if (exist('BoutSpikeChanNo', 'var'))
                [TempData, Fs] = SSAReadOKrankData(DataDirectory, DataDirectory, Neural_INR.BoutDetails(i).SongFile, BoutSpikeChanNo(i));
            else
                [TempData, Fs] = SSAReadOKrankData(DataDirectory, DataDirectory, Neural_INR.BoutDetails(i).SongFile, SpikeChanNo);
            end
        else
           if (strfind(FileType, 'obs'))
                channel_string = strcat('obs',num2str(5 - SpikeChanNo),'r');
                [TempData, Fs] = SSASoundIn([DataDirectory, '/'], DataDirectory, Neural_INR.BoutDetails(i).SongFile, channel_string);
                % Convert to V - 5V on the data acquisition is 32768
                TempData = TempData * 5/32768;
            end
        end
        [b, a] = butter(4, [600/16000 6000/16000]);
        TempData = filtfilt(b, a, TempData);         
        TempData = TempData(:);
        if (exist('BoutGain', 'var'))
            TempData = TempData*BoutGain(i);
        end
        SpikeRawData = [SpikeRawData; TempData];

        SpikeIndices = find((Neural_INR.BoutDetails(i).SpikeTimes >= (Neural_INR.BoutDetails(i).onsets(1) - 0.05)) & (Neural_INR.BoutDetails(i).SpikeTimes >= (Neural_INR.BoutDetails(i).offsets(1))));
        SpikeIndices = round(Neural_INR.BoutDetails(i).SpikeTimes * Fs);
        SpikeIndices = SpikeIndices(:);
        SpikeWaveformIndices = repmat(SpikeIndices, 1, length(WFIndices));
        SpikeWaveformIndices = SpikeWaveformIndices + repmat(WFIndices, size(SpikeWaveformIndices,1), 1);
        if (size(SpikeWaveformIndices,1) == 1)
            Waveforms = [Waveforms; SpikeRawData(SpikeWaveformIndices)'];
        else
            Waveforms = [Waveforms; SpikeRawData(SpikeWaveformIndices)];
        end
    end
end

if (InvertWaveform == 1)
    Waveforms = -Waveforms;
end

AlignedWaveforms = AlignSpikeWaveformsByMax(Waveforms, 3, Threshold, PreSamples);

MeanSpikeWaveform = mean(AlignedWaveforms);
MeanSubractedWaveforms = AlignedWaveforms - repmat(MeanSpikeWaveform, size(AlignedWaveforms,1), 1);

SignalToNoise.SongSignalToNoise = (max(MeanSpikeWaveform) - min(MeanSpikeWaveform))/std(MeanSubractedWaveforms(:)*2);
SignalToNoise.OneSideSongSignalToNoise = [max(MeanSpikeWaveform)/std(MeanSubractedWaveforms(:)) min(MeanSpikeWaveform)/std(MeanSubractedWaveforms(:))];

if (~iscell(FileList))
    TempFileList = FileList;
    clear FileList;
    for i = 1:size(TempFileList, 1)
        FileList{i} = TempFileList(i,:);
    end
end
 
AllWaveforms = [];
for i = 1:length(FileList),
    if (strfind(FileType, 'okrank'))
        if (exist('BoutSpikeChanNo', 'var'))
            [TempData, Fs] = SSAReadOKrankData(DataDirectory, DataDirectory, FileList{i}, FileSpikeChanNo(i));
        else
            [TempData, Fs] = SSAReadOKrankData(DataDirectory, DataDirectory, FileList{i}, SpikeChanNo);
        end
    else
       if (strfind(FileType, 'obs'))
            channel_string = strcat('obs',num2str(5 - SpikeChanNo),'r');
            [TempData, Fs] = SSASoundIn([DataDirectory, '/'], DataDirectory, FileList{i}, channel_string);
            % Convert to V - 5V on the data acquisition is 32768
            TempData = TempData * 5/32768;
        end
    end
    [b, a] = butter(4, [300/16000 8000/16000]);
    TempData = filtfilt(b, a, TempData);
    if (exist('FileBoutGain', 'var'))
        TempData = TempData*FileBoutGain(i);
    end
    
    TempSpikeTimes = load([SpikeDir, '/Chan', num2str(SpikeChanNo), '_', FileList{i}, '.spk']);
    SpikeIndices = [];
    for i = 1:length(ClusterNos),
        SpikeIndices = [SpikeIndices; find(TempSpikeTimes(:,1) == ClusterNos(i))];
    end
    SpikeIndices = round(TempSpikeTimes(SpikeIndices,2) * Fs);
    SpikeIndices = SpikeIndices(:);
    
    RemovalIndices = find(SpikeIndices <= abs(WFIndices(1)));
    SpikeIndices(RemovalIndices) = [];
    
    RemovalIndices = find(SpikeIndices >= (length(TempData) - WFIndices(end)));
    SpikeIndices(RemovalIndices) = [];
    
    SpikeWaveformIndices = repmat(SpikeIndices, 1, length(WFIndices));
    SpikeWaveformIndices = SpikeWaveformIndices + repmat(WFIndices, size(SpikeWaveformIndices,1), 1);

    if (size(SpikeWaveformIndices,1) == 1)
       AllWaveforms = [AllWaveforms; TempData(SpikeWaveformIndices)'];
    else
       AllWaveforms = [AllWaveforms; TempData(SpikeWaveformIndices)];
    end
end

if (InvertWaveform == 1)
    AllWaveforms = -AllWaveforms;
end

AllAlignedWaveforms = AlignSpikeWaveformsByMax(AllWaveforms, 3, Threshold, PreSamples);

MeanAllSpikeWaveform = mean(AllAlignedWaveforms);
AllResidualWaveforms = AllAlignedWaveforms - repmat(MeanAllSpikeWaveform, size(AllAlignedWaveforms,1), 1);

SignalToNoise.AllSpikeSignalToNoise = (max(MeanAllSpikeWaveform) - min(MeanAllSpikeWaveform))/(2*std(AllResidualWaveforms(:)));
SignalToNoise.OneSideAllSpikeSignalToNoise = [max(MeanAllSpikeWaveform)/std(AllResidualWaveforms(:)) min(MeanAllSpikeWaveform)/std(AllResidualWaveforms(:))];
figure;
subplot(1,2,1);
plot(AlignedWaveforms', 'Color', 'b');
hold on;
plot(MeanSpikeWaveform, 'r');
axis tight;
temp = axis;
text(120, temp(4)/2, ['SNR = ', num2str(SignalToNoise.SongSignalToNoise)]);
title('Song waveforms');

subplot(1,2,2);
plot(AllAlignedWaveforms', 'Color', 'b');
hold on;
plot(MeanAllSpikeWaveform, 'r');
axis tight;
temp = axis;
text(120, temp(4)/2, ['SNR = ', num2str(SignalToNoise.AllSpikeSignalToNoise)]);
title('All waveforms');

%saveas(gcf, ['/home/raghav/HVCPlots/', FileList{1}, '.SpikeSNR.tiff'], 'tiff');

disp('Finished');