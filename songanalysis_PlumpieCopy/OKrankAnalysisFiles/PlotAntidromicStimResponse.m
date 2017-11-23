function [] = PlotAntidromicStimResponse(DirectoryName, FileName, SpikeChanNo, StimChanNo, TitleString, PreTime, PostTime, Threshold, EvokedSpike_ApproxTime, MinorMax, varargin)

if (nargin > 10)
    StimtoPlot = varargin{1};
end

[SpikeChan, Fs] = ReadOKrankData(DirectoryName, FileName, SpikeChanNo);
[b, a] = butter(3, [600/16000 6000/16000]);
SpikeChan = filtfilt(b, a, SpikeChan);
SpikeChan = SpikeChan * 100;
[StimChan, Fs] = ReadOKrankData(DirectoryName, FileName, StimChanNo);

h = [1 -1];
Temp = StimChan >= Threshold;
temp = zeros(size(Temp, 1), size(Temp, 2));
temp(find(Temp > 0)) = 1;
trans = conv(h, temp);
onsets = find(trans > 0);

PreNSamples = round(PreTime * Fs/1000);
PostNSamples = round(PostTime * Fs/1000);

Time = (1:1:(PreNSamples + PostNSamples + 1))/Fs - PreNSamples/Fs;
Time = Time * 1000;
MoreThan_x_ms = find(Time >= EvokedSpike_ApproxTime);

Index = 1;
for i = 1:2:length(onsets),
    if ((onsets(i) - PreNSamples) < 1)
        continue;
    end
    if ((onsets(i) + PostNSamples) > length(SpikeChan))
        continue;
    end
    EvokedWaveform(Index,:) = SpikeChan((onsets(i) - PreNSamples):(onsets(i) + PostNSamples));
    StimChanWaveform(Index,:) = StimChan((onsets(i) - PreNSamples):(onsets(i) + PostNSamples));
    if (strfind(MinorMax, 'min'))
        [MaxVal, MaxValInd] = min(EvokedWaveform(Index,MoreThan_x_ms));
    else
        [MaxVal, MaxValInd] = max(EvokedWaveform(Index,MoreThan_x_ms));
    end
    AntidromicSpikeTime(Index) = Time(MaxValInd - 1 + MoreThan_x_ms(1));
    AntidromicSpikeIndex(Index) = MaxValInd - 1 + MoreThan_x_ms(1);
    Index = Index + 1;
end

if (exist('StimtoPlot', 'var'))
    disp(['No of trials is ', num2str(length(StimtoPlot))]);
    disp(['Mean spike time after antidromic stimulation is ', num2str(mean(AntidromicSpikeTime(StimtoPlot))), ' ms and std is ', num2str(std(AntidromicSpikeTime(StimtoPlot))), ' ms']);
    disp(['Mean spike time after antidromic stimulation is ', num2str(mean(AntidromicSpikeTime(StimtoPlot))), ' ms and 10-90th percentile is ', num2str(prctile(AntidromicSpikeTime(StimtoPlot), 90) - prctile(AntidromicSpikeTime(StimtoPlot),10)), ' ms']);
else
    disp(['No of trials is ', num2str(Index - 1)]);
    disp(['Mean spike time after antidromic stimulation is ', num2str(mean(AntidromicSpikeTime)), ' ms and std is ', num2str(std(AntidromicSpikeTime)), ' ms']);
    disp(['Mean spike time after antidromic stimulation is ', num2str(mean(AntidromicSpikeTime)), ' ms and std by percentile is ', num2str((prctile(AntidromicSpikeTime, 84) - prctile(AntidromicSpikeTime,16))/2), ' ms']);
end
for i = 1:Index-1,
    if (mod(i,16) == 1)
        figure;
        SubPlotIndex = 1;
    end
    subplot(4,4,SubPlotIndex);
    plot(EvokedWaveform(i,:));
    hold on;
    plot(StimChanWaveform(i,:)*100, 'r');
    SubPlotIndex = SubPlotIndex + 1;
end

figure
if (exist('StimtoPlot', 'var'))
    plot(Time, EvokedWaveform(StimtoPlot,:)');
    hold on;
    for i = 1:length(StimtoPlot),
        plot(Time(AntidromicSpikeIndex(StimtoPlot(i))), EvokedWaveform(i,AntidromicSpikeIndex(StimtoPlot(i))), 'rs');
    end
else
    plot(Time, EvokedWaveform', 'k');
    hold on;
    for i = 1:length(AntidromicSpikeIndex),
        plot(Time(AntidromicSpikeIndex(i)), EvokedWaveform(i,AntidromicSpikeIndex(i)), 'rs');
    end
end
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gcf, 'Color', 'w');
xlabel('Time (ms)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Amplitude (uV)', 'FontSize', 12, 'FontWeight', 'bold');
axis tight
title([FileName, ' ', TitleString], 'FontSize', 12, 'FontWeight', 'bold');