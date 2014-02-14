function [DNoteInfo, UNoteInfo] = PlotSegmentedSAPFeatures(DirectoryName, FileType)

cd(DirectoryName);
!mkdir RecFiles

NoteFiles = dir('*_dir*.');
NoteInfo = [];
AvgPower = [];
for i = 1:length(NoteFiles),
    disp(NoteFiles(i).name);
    try
        SAPFeatures = CalculateSAPFeatures(NoteFiles(i).name(1:(end-8)), 44100);
    catch
        disp(['Could not calculate SAP Features']);
        continue;
    end
    
    if (strfind(FileType, 'okrank'))
        [RawData, Fs] = ReadOKrankData(handles.MakeTemplatesData.DataDirectory, handles.MakeTemplatesData.SongFileName, 0);
else
    if (strfind(handles.MakeTemplatesData.FileType, 'obs'))
        [RawData, Fs] = soundin_copy(handles.MakeTemplatesData.DataDirectory, handles.MakeTemplatesData.SongFileName, handles.MakeTemplatesData.FileType);   
        RawData = RawData/32768;
    else
        if (strfind(handles.MakeTemplatesData.FileType, 'wav'))
            [RawData, Fs] = wavread(handles.MakeTemplatesData.SongFileName);    
        end
    end
end

Time = (1:1:length(RawData))/Fs;

% Now using an 8 pole butterworth bandpass filter as default.
[b,a]=butter(8,[300*2/Fs, 10000*2/Fs]);

FiltSong=filtfilt(b, a, RawData);
  
if length(RawData) ~= length(FiltSong) 
  disp(['warning! bandpass: input and output file lengths do not match!']);
end

RawData = FiltSong;
clear FiltSong;

[Spect, Freq, SpectTime, Power] = spectrogram(RawData, handles.MakeTemplatesData.FFTWindowLength, handles.MakeTemplatesData.FFTWindowOverlap, handles.MakeTemplatesData.FFTWindowLength, Fs);
FreqRows = find((Freq >= 1700) & (Freq <= 7100));

Spect = Spect(FreqRows,:);
Power = Power(FreqRows,:);

axes(handles.AmplitudePlot);
hold off;
%plot(Time, RawData);
plot(SpectTime, 10*log10(sum(Power)), 'r');
hold on;

axis tight;

set(handles.AmplitudePlot, 'FontSize', 12, 'FontWeight', 'bold');
%xlabel('Time (sec)', 'FontSize', 12, 'FontWeight', 'bold');
%ylabel('Amplitude', 'FontSize', 12, 'FontWeight', 'bold');

handles.MakeTemplatesData.AmpPlotAxis = axis;
handles.MakeTemplatesData.AmpPlotAxis(1:2) = [Time(1) Time(end)];

SpectFs = 1/(SpectTime(2) - SpectTime(1));

[onsets, offsets] = segment(10*log10(sum(Power)), SpectFs, 5, 10, handles.MakeTemplatesData.Threshold);

onsets = onsets/1000;
offsets = offsets/1000;

    %Notes = load(NoteFiles(i).name);
    Notes.onsets = Notes.onsets/1000;
    Notes.offsets = Notes.offsets/1000;
    for j = 1:length(Notes.onsets),
        NoteInfo((end+1),1) = (Notes.offsets(j) - Notes.onsets(j)) * 1000;
        StartIndex = find(SAPFeatures.Time < Notes.onsets(j), 1, 'last');
        EndIndex = find(SAPFeatures.Time < Notes.offsets(j), 1, 'last');
        NoteInfo(end, 2) = mean(SAPFeatures.AM(StartIndex:EndIndex));
        NoteInfo(end, 3) = mean(SAPFeatures.Entropy(StartIndex:EndIndex));
        NoteInfo(end, 4) = mean(SAPFeatures.Freq(StartIndex:EndIndex));
        NoteInfo(end, 5) = mean(SAPFeatures.Pitch(StartIndex:EndIndex));
        NoteInfo(end, 6) = mean(SAPFeatures.FM(StartIndex:EndIndex));
        NoteInfo(end, 7) = mean(SAPFeatures.Amplitude(StartIndex:EndIndex));
        NoteInfo(end, 8) = mean(SAPFeatures.PitchGoodness(StartIndex:EndIndex));
    end
end
DNoteInfo = NoteInfo;

NoteFiles = dir('*_undir*.not.mat');
NoteInfo = [];
for i = 1:length(NoteFiles),
    disp(NoteFiles(i).name);
    try
        SAPFeatures = CalculateSAPFeatures(NoteFiles(i).name(1:(end-8)), 44100);
    catch
        disp(['Could not calculate SAP Features']);
        continue;
    end
    
    Notes = load(NoteFiles(i).name);
    Notes.onsets = Notes.onsets/1000;
    Notes.offsets = Notes.offsets/1000;
    for j = 1:length(Notes.onsets),
        NoteInfo((end+1),1) = (Notes.offsets(j) - Notes.onsets(j)) * 1000;
        StartIndex = find(SAPFeatures.Time < Notes.onsets(j), 1, 'last');
        EndIndex = find(SAPFeatures.Time < Notes.offsets(j), 1, 'last');
        NoteInfo(end, 2) = mean(SAPFeatures.AM(StartIndex:EndIndex));
        NoteInfo(end, 3) = mean(SAPFeatures.Entropy(StartIndex:EndIndex));
        NoteInfo(end, 4) = mean(SAPFeatures.Freq(StartIndex:EndIndex));
        NoteInfo(end, 5) = mean(SAPFeatures.Pitch(StartIndex:EndIndex));
        NoteInfo(end, 6) = mean(SAPFeatures.FM(StartIndex:EndIndex));
        NoteInfo(end, 7) = mean(SAPFeatures.Amplitude(StartIndex:EndIndex));
        NoteInfo(end, 8) = mean(SAPFeatures.PitchGoodness(StartIndex:EndIndex));
    end
end
UNoteInfo = NoteInfo;
