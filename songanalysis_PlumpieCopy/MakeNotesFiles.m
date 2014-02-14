function [] = MakeNotesFiles(DirectoryName, FileName, FileType, Threshold)

cd(DirectoryName);

if (strfind(FileType, 'obs'))
    [RawData, Fs] = soundin_copy(DirectoryName, FileName, 'obs0r');
else
    if (strfind(FileType, 'okrank'))
        [RawData,Fs] = ReadOKrankData(DirectoryName, FileName, 0);
    end
end

Time = (1:1:length(RawData))/Fs;

figure
plot(Time, RawData);
hold on;
[b, a] = butter(8, [300*2/Fs, 10000*2/Fs]);

FiltSong = filtfilt(b, a, RawData);

[S, F, T, P] = spectrogram(RawData, 128, 124, 128, Fs, 'yaxis');
FreqRows = find((F >= 1500) & (F <= 7300));
AmpEnvelope = P(FreqRows, :);
AmpEnvelope = AmpEnvelope - mean(mean(AmpEnvelope));
AmpEnvelope = AmpEnvelope/(sqrt((sum(sum(AmpEnvelope.*AmpEnvelope)))/(size(AmpEnvelope,1) * size(AmpEnvelope,2))));
AmpEnvelope = sum(AmpEnvelope);
AmpEnvelope = AmpEnvelope/max(AmpEnvelope);
plot(T, AmpEnvelope, 'r');
axis tight;
temp = axis;

NewFs = 1/(T(2) - T(1));
[onsets, offsets] = segment(AmpEnvelope, NewFs, 5, 10, Threshold);
for i = 1:length(onsets),
    plot([onsets(i)/1000 onsets(i)/1000 offsets(i)/1000 offsets(i)/1000], [0 0.9*temp(4) 0.9*temp(4) 0], 'm');
end
min_int = 5;
min_dur = 10;

for i = 1:length(onsets),
    labels(i) = '0';
end
sm_win = 2;
Window = 128;
NFFT = 128;
NOverlap = 124;

%save([FileName, '.not.mat'], 'Fs', 'labels', 'min_int', 'min_dur', 'offsets', 'onsets', 'sm_win', 'Threshold', 'Window', 'NFFT', 'NOverlap');
