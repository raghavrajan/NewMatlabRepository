function [] = MatchTemplates(DirectoryName, SongFileName, TemplateFileName, FileType, TemplateThreshold, FinalThreshold, PlotOption)

FFTWindowLength = 256;
FFTWindowOverlap = 128;

load(TemplateFileName);

if (ispc)
    if (DirectoryName(end) ~= '\')
        DirectoryName(end + 1) = '\';
    end
else
    if (DirectoryName(end) ~= '/')
        DirectoryName(end + 1) = '/';
    end 
end
   
if (strfind(FileType, 'okrank'))
        [RawData, Fs] = ReadOKrankData(DirectoryName, SongFileName, 0);
else
    if (strfind(FileType, 'obs'))
        [RawData, Fs] = soundin_copy(DirectoryName, SongFileName, FileType);  
        RawData = RawData/32768;
    else
        if (strfind(FileType, 'wav'))
            [RawData, Fs] = wavread(SongFileName);    
        end
    end
end

% Now using an 8 pole butterworth bandpass filter as default.
[b,a]=butter(8,[300*2/Fs, 10000*2/Fs]);

FiltSong=filtfilt(b, a, RawData);
  
if length(RawData) ~= length(FiltSong) 
  disp(['warning! bandpass: input and output file lengths do not match!']);
end

RawData = FiltSong;
clear FiltSong;

[Spect, Freq, SpectTime, Power] = spectrogram(RawData, FFTWindowLength, FFTWindowOverlap, FFTWindowLength, Fs);
FreqRows = find((Freq >= 1700) & (Freq <= 7100));

Spect = Spect(FreqRows, :);
Power = Power(FreqRows, :);

if (strfind(PlotOption, 'yes'))
    figure;
    SpectrogramPlot = subplot(1,1,1);
    plot_motif_spectrogram(RawData, Fs, gcf, SpectrogramPlot);
    hold on;
end

NewFs = 1/(SpectTime(2) - SpectTime(1));
[onsets, offsets] = segment(sum(abs(Spect)), NewFs, 5, 10, TemplateThreshold);
if (isempty(onsets))
    disp('No syllables in this file');
    return;
end
onsets = onsets/1000;
offsets = offsets/1000;

Colours = ['g' 'b' 'c' 'm' 'y' 'k'];

for i = 1:length(Templates),
    %clear Correlation2;
    for j = 1:length(Templates(i).STemplate),
        Index = 1;
        while(Index <= (size(Spect,2) - size(Templates(i).STemplate{j},2))),
            %if (length(find(onsets == SpectTime(Index))) > 0)
            %    disp('Reached a syllable');
            %end
            Temp = Spect(:,(Index:(Index + size(Templates(i).STemplate{j},2) - 1)));
            Temp = ScaleSpect(abs(Temp));
            MeanTemp = sum(sum(Temp))/((size(Temp, 1) * size(Temp, 2)));
            STDTemp = sqrt((sum(sum((Temp - MeanTemp).*(Temp - MeanTemp))))/((size(Temp, 1) * size(Temp, 2)) - 1));
            Temp = (Temp - MeanTemp)/STDTemp;
            STemplateMatch(i).Value{j}(Index) = 1/(sum(sum(abs(Templates(i).STemplate{j} - (Temp)))));
            
            Index = Index + 1;
        end
        STemplateMatch(i).Value{j} = (STemplateMatch(i).Value{j} - mean(STemplateMatch(i).Value{j}))/(std(STemplateMatch(i).Value{j}));
    end
end


for i = 1:length(STemplateMatch),
    MaxMatch(i) = 0;
    for j = 1:length(STemplateMatch(i).Value),
        MaxMatch(i) = MaxMatch(i) + max(STemplateMatch(i).Value{j});
    end
    MaxMatch(i) = MaxMatch(i)/j;
end

[MaxValue, MaxIndex] = max(MaxMatch);
for i = 1:length(STemplateMatch),
    for j = 1:length(STemplateMatch(i).Value),
        [Peaks{i, j}, Locations{i, j}] = findpeaks(STemplateMatch(i).Value{j},'minpeakheight', MaxMatch(i)/2);
    end
end

[Spect, Freq, SpectTime1, Power] = spectrogram(RawData, 128, 124, 128, Fs);
FreqRows = find((Freq >= 1700) & (Freq <= 7100));

Spect = Spect(FreqRows, :);
Power = Power(FreqRows, :);

NewFs = 1/(SpectTime1(2) - SpectTime1(1));
[NewOnsets, NewOffsets] = segment(sum(abs(Spect)), NewFs, 5, 10, FinalThreshold);
NewOnsets = NewOnsets/1000;
NewOffsets = NewOffsets/1000;

for i = 1:size(Locations,1),
    MinIndex = 1000;
    for j = 1:size(Locations,2),
        if (length(Locations{i,j}) ~= 0)
            if (length(Locations{i,j}) < MinIndex)
                MinIndex = length(Locations{i,j});
                ColNo = j;
            end
        end
    end
    Syllables{i}.Location(ColNo,:) = [Locations{i,ColNo}];
    for j = 1:size(Locations, 2),
        if (j == ColNo),
            continue;
        end
        if (length(Locations{i,j}) == 0)
            continue;
        end
        for k = 1:MinIndex,
            TempIndex = find((Locations{i,j} > (Syllables{i}.Location(ColNo,k) - 5)) & (Locations{i,j} < (Syllables{i}.Location(ColNo,k) + 5)));
            if (isempty(TempIndex))
                Syllables{i}.Location(j,k) = Syllables{i}.Location(ColNo,k);
            else
                Syllables{i}.Location(j,k) = Locations{i,j}(TempIndex(1));
            end
        end
    end
end

Index = 1;
clear onsets;
clear offsets;
clear labels;

for i = 1:length(Syllables),
    Syllables{i}.Location = round(mean(Syllables{i}.Location));

    OnsetTimes = SpectTime(Syllables{i}.Location);
    
    OnsetIndex = [];
    for j = 1:length(OnsetTimes),
        TempOnset1 = find(NewOnsets < OnsetTimes(j), 1, 'last');
        TempOnset2 = find(NewOnsets > OnsetTimes(j), 1, 'first');
        
        if ((~isempty(TempOnset1)) && (~isempty(TempOnset2)))
            if ((OnsetTimes(j) - NewOnsets(TempOnset1)) > (NewOnsets(TempOnset2) - OnsetTimes(j)))
                onsets(Index) = NewOnsets(TempOnset2);
                OnsetIndex = TempOnset2;
            else
                onsets(Index) = NewOnsets(TempOnset1);
                OnsetIndex = TempOnset1;
            end
        else
            if (~isempty(TempOnset1))
                onsets(Index) = NewOnsets(TempOnset1);
                OnsetIndex = TempOnset1;
            else
                if (~isempty(TempOnset2))
                    onsets(Index) = NewOnsets(TempOnset2);
                    OnsetIndex = TempOnset2;
                else
                    onsets(Index) = -5;
                    offsets(Index) = -5;
                    labels(Index) = '0';
                end
            end
        end
        if (onsets(Index) > 0)
            OffsetIndex = OnsetIndex;
            TempDuration = NewOffsets(OffsetIndex) - onsets(Index);
            while (~((TempDuration > (0.7 * Templates(i).Duration)) && (TempDuration < (1.3 * Templates(i).Duration))))
                if (TempDuration < (0.75 * Templates(i).Duration))
                    OffsetIndex = OffsetIndex + 1;
                    if (OffsetIndex > length(NewOffsets))
                        OffsetIndex = OffsetIndex - 1;
                        break;
                    end
                    if (OffsetIndex > (OnsetIndex + 3))
                        break;
                    end
                else
                    if (OffsetIndex ~= OnsetIndex)
                        OffsetIndex = OffsetIndex -1;
                    end
                    break;
                end
                TempDuration = NewOffsets(OffsetIndex) - onsets(Index);
            end
            offsets(Index) = NewOffsets(OffsetIndex);
            labels(Index) = Templates(i).Label;
        end
        Index = Index + 1;
    end
end

if (~isempty(onsets))
    Indices = find(onsets < 0);
    onsets(Indices) = [];
    offsets(Indices) = [];
    labels(Indices) = [];

    [SortedValues, SortedIndices] = sort(onsets);
    onsets = onsets(SortedIndices);
    offsets = offsets(SortedIndices);
    labels = labels(SortedIndices);

    [UniqueOnsets, NewIndices, OrigIndices] = unique(onsets);
    if (length(UniqueOnsets) ~= length(onsets))
        for i = 1:OrigIndices(end),
            Indices = find(OrigIndices == i);
            if (length(i) > 1)
                SyllableTime = onsets(Indices(1));
                DuplicateLabels = labels(Indices);
                clear TemplateIndices;
                TemplateIndices = [];
                for j = 1:length(DuplicateLabels),
                    TemplateIndices = [TemplateIndices find([Templates.Label] == DuplicateLabels(j))];
                end
                TimeIndices = find((SpectTime > (SyllableTime - 0.005)) & (SpectTime < (SyllableTime + 0.005)));
                for j = 1:length(STemplateMatch),
                    if (isempty(find(TimeIndices == j)))
                        continue;
                    end
                    for k = 1:length(STemplateMatch(j).Value),
                       Values{j} 
    
    if (strfind(PlotOption, 'yes'))
        for i = 1:length(onsets),
            plot([onsets(i) onsets(i) offsets(i) offsets(i)], [0 5000 5000 0], 'r');
            PlotLabels(i) = text([onsets(i)], 5100, labels(i));
            set(PlotLabels(i), 'FontSize', 16, 'FontWeight', 'bold');
        end
    end

    onsets = onsets * 1000;
    offsets = offsets * 1000;
    sm_win = 2;
    threshold = 0.0005;
    min_dur = 10;
    min_int = 5;
    save([SongFileName, '.not.mat'], 'Fs', 'labels', 'onsets', 'offsets', 'min_dur', 'min_int', 'offsets', 'onsets', 'sm_win', 'threshold');
    disp(labels);
else
    disp('No syllable matches in this file');
end

