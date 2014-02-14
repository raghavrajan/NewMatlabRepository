function [] = AutoSegmentAndLabel(DirectoryName, FileName, FileType, SongChanNo, TemplateFile, Threshold)

Templates = load(TemplateFile);

FFTWinSize = 0.005; % in sec
FFTWinOverlap = 0.95;

if (strfind(FileType, 'okrank'))
    [Data, Fs] = ReadOKrankData(DirectoryName, FileName, SongChanNo);
end
Time = (1:1:length(Data))/Fs;

WinSize = round(FFTWinSize * Fs);
WinOverlap = round(FFTWinOverlap * WinSize);

[S, F, T, P] = spectrogram(Data, hamming(WinSize), WinOverlap, WinSize, Fs);

Amp = 10*log10(sum(P));
NoteTimes = Amp > Threshold;
h = [1 -1];
TempNoteTimes = zeros(size(NoteTimes));
TempNoteTimes(NoteTimes > 0) = 1;
NoteTimes = conv(TempNoteTimes, h, 'same');
clear TempNoteTimes;

OnsetIndices = find(NoteTimes > 0);
OffsetIndices = find(NoteTimes < 0);

if (length(OnsetIndices) > length(OffsetIndices))
    OnsetIndices(end) = [];
else
    if (length(OffsetIndices) > length(OnsetIndices))
        OffsetIndices(end) = [];
    end
end
Onsets = T(OnsetIndices);
Offsets = T(OffsetIndices);

Durations = Offsets - Onsets;
ShortNotes = Durations < 0.015;
Onsets(ShortNotes) = [];
Offsets(ShortNotes) = [];
OnsetIndices(ShortNotes) = [];
OffsetIndices(ShortNotes) = [];
Durations = Offsets - Onsets;

figure(1);
plot(T, Amp);
hold on;
AmpFs = 1/(T(2) - T(1));
SpectWin = round(0.002 * AmpFs);

for i = 1:length(Onsets),
    OnsetIndices(i) = find(Time <= Onsets(i), 1, 'last');
    OffsetIndices(i) = find(Time <= Offsets(i), 1, 'last');
end

SpectWin = round(0.005 * Fs);
SpectOverlap = round(0.9 * SpectWin);
FreqVect = [300:75:10000];

for i = 1:length(Onsets),
     figure(1);
     plot([Onsets(i) Onsets(i) Offsets(i) Offsets(i)], [-70 0 0 -70], 'r');
     [Spec{i} Freq{i} SpecTime{i} Power{i}] = spectrogram(Data(OnsetIndices(i):OffsetIndices(i)), [], 0, FreqVect, Fs, 'yaxis');
end

Colours = ['rgbcmyk'];
for i = 1:length(Templates.Labels),
    Template = Templates.PSD{i};
    for j = 1:length(Onsets),
        Template1 = Power{j};
        MeanTemp1 = sum(sum(Template1))/(size(Template1, 1) * size(Template1, 2));
        STDTemp1 = sqrt((sum(sum((Template1 - MeanTemp1).*(Template1 - MeanTemp1))))/((size(Template1, 1) * size(Template1,2)) - 1));
        Template1 = (Template1 - MeanTemp1)/STDTemp1;
        Diff(i,j) = 1/(sum(sum(abs(Template1 - Template))));
        if (Templates.Durations{i} >= Durations(j))
            Diff(i,j) = Diff(i,j) * Durations(j)/Templates.Durations{i};
        else
            Diff(i,j) = Diff(i,j) * Templates.Durations{i}/Durations(j);
        end
    end
    plot(Onsets, ((Diff(i,:)/max(Diff(i,:)) * 20)-60), [Colours(mod(i-1,7) + 1), 's-']);
end

PlotSpectrogram(DirectoryName, FileName, FileType);
hold on;
for i = 1:length(Onsets),
    [MaxValue, MaxIndex] = max(Diff(:,i));
    text((Offsets(i) + Onsets(i))/2, 8000, Templates.Labels(MaxIndex), 'FontSize', 12, 'FontWeight', 'bold');
end
disp('Finished');
