function [] = PlotAntidromicStimResponse(DirectoryName, FileName, SpikeChanNo, StimChanNo, TitleString, PreTime, PostTime, Threshold, varargin)

if (nargin > 8)
    StimtoPlot = varargin{1};
end

[SpikeChan, Fs] = ReadOKrankData(DirectoryName, FileName, SpikeChanNo);
[b, a] = butter(3, [600/16000 6000/16000]);
SpikeChan = filtfilt(b, a, SpikeChan);
SpikeChan = SpikeChan * 100;
[StimChan, Fs] = ReadOKrankData(DirectoryName, FileName, StimChanNo);

h = [1 -1];
Temp = StimChan > Threshold;
temp = zeros(size(Temp, 1), size(Temp, 2));
temp(find(Temp > 0)) = 1;
trans = conv(h, temp);
onsets = find(trans > 0);

Index = 1;
PreNSamples = round(PreTime * Fs/1000);
PostNSamples = round(PostTime * Fs/1000);
for i = 1:2:length(onsets),
    if ((onsets(i) - PreNSamples) < 1)
        continue;
    end
    if ((onsets(i) + PostNSamples) > length(SpikeChan))
        continue;
    end
    EvokedWaveform(Index,:) = SpikeChan((onsets(i) - PreNSamples):(onsets(i) + PostNSamples));
    Index = Index + 1;
end

disp(['No of trials is ', num2str(Index - 1)]);
for i = 1:Index-1,
    if (mod(i,16) == 1)
        figure;
        SubPlotIndex = 1;
    end
    subplot(4,4,SubPlotIndex);
    plot(EvokedWaveform(i,:));
    SubPlotIndex = SubPlotIndex + 1;
end

figure
Time = (1:1:size(EvokedWaveform,2))/Fs - PreNSamples/Fs;
Time = Time * 1000;
if (exist('StimtoPlot', 'var'))
    plot(Time, EvokedWaveform(StimtoPlot,:)');
else
    plot(Time, EvokedWaveform', 'k');
end
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gcf, 'Color', 'w');
xlabel('Time (ms)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Amplitude (uV)', 'FontSize', 12, 'FontWeight', 'bold');
axis tight
title([FileName, ' ', TitleString], 'FontSize', 12, 'FontWeight', 'bold');