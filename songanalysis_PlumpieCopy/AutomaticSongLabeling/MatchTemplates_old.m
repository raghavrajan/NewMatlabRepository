function [] = MatchTemplates(DirectoryName, SongFileName, TemplateFileName, FileType)

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
        RawData = RawData/max(RawData);
        RawData = RawData - mean(RawData);
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

Freq  = Fs/2*linspace(0,1,512/2+1);
FreqRows = find((Freq >= 1700) & (Freq <= 7100));
[Y, F, T, P] = spectrogram(RawData, 512, 256, 512, Fs, 'yaxis');
Time = (1:1:length(RawData))/Fs;

nfft=round(Fs*2/1000);
nfft = 2^nextpow2(nfft);
spect_win = hanning(nfft);
noverlap = round(0.5*length(spect_win)); %number of overlapping points       
%now calculate spectrogram
[spect, freq, time, power] = spectrogram(RawData,spect_win,noverlap,nfft,Fs,'yaxis');
T = time;
FreqRows = find((freq >= 1500) & (freq <= 7300));
%FreqRows = find((F >= 1500) & (F <= 7300));
idx_spect=scale_spect(spect);  %calculate index array for spectrogram
P = idx_spect;
%P = 20*log10(Y);

% subplot(2,1,2);
% plot(Time,RawData);
% axis tight;
figure;
SpectrogramPlot = subplot(1,1,1);
plot_motif_spectrogram(RawData, Fs, gcf, SpectrogramPlot);
hold on;
plot(T, sum(power)/max(sum(power)) * 5000, 'r');
AmplitudeEnvelope = abs(sum(power)/max(sum(power)));
NewFs = 1/(T(2) - T(1));
[onsets, offsets] = segment(AmplitudeEnvelope, NewFs, 5, 25, 0.0005);
onsets = onsets/1000;
offsets = offsets/1000;

for i = 1:length(onsets),
    plot([onsets(i) onsets(i) offsets(i) offsets(i)], [0 5000 5000 0], 'r');
end

Colours = ['g' 'b' 'c' 'm' 'y' 'k'];

for i = 1:length(Templates),
%    if (mod(i,3) == 1)
%        figure;
%    end
    Index = 1;
    clear Correlation Correlation2 Correlation3 TemplateMatchValue TemplateMatchValue2 TemplateMatchValue3;
    while(Index <= (size(P,2) - size(Templates(i).Template,2))),
        Temp = P(:,(Index:(Index + size(Templates(i).Template,2) - 1)));
        Temp = Temp - mean(mean(Temp));
        Temp = Temp/(sqrt((sum(sum(Temp.*Temp)))/(size(Temp,1) * size(Temp,2))));
        
        Temp2 = sum(Temp(FreqRows,:));
        TemplateMatchValue(Index) = 1/(sum(sum(abs(Templates(i).Template - Temp))));
        TemplateMatchValue2(Index) = 1/(sum(sum(abs(Templates(i).Template(FreqRows,:) - Temp(FreqRows,:)))));
        Correlation(Index) = 1/(sum(abs(sum(Templates(i).Template) - sum(Temp))));
        Correlation2(Index) = 1/(sum(abs(sum(Templates(i).Template(FreqRows,:)) - sum(Temp(FreqRows,:)))));
        Index = Index + 1;
    end

    Index = 1;
    while(Index <= (size(P,2) - size(Templates(i).Template2,2))),
        Temp = P(FreqRows,(Index:(Index + size(Templates(i).Template2,2) - 1)));
        Temp = Temp - mean(mean(Temp));
        Temp = Temp/(sqrt((sum(sum(Temp.*Temp)))/(size(Temp,1) * size(Temp,2))));
        TemplateMatchValue3(Index) = 1/(sum(sum(abs(Templates(i).Template2 - Temp))));
        Correlation3(Index) = 1/(sum(abs(sum(Templates(i).Template2) - sum(Temp))));
        Index = Index + 1;
    end
    TemplateMatch{i}.Value = TemplateMatchValue;
    [Test, TestIndices] = sort(Correlation2);
    Indices = find(Correlation2 > Test(round(0.985 * length(Test))));
    [pks, locs] = findpeaks(Correlation2(Indices));
    %SpectrogramPlot = subplot(1,1,1);
    %MedianMotif.FileName = SongFileName;
    %MedianMotif.Length = Time(end);
    %MedianMotif.StartTime = Time(1);

    %PlotMotifSpectrogram(DirectoryName, FileType, MedianMotif, gcf, subplot(length(Templates),1,i));

%    plot_motif_spectrogram(RawData, Fs, gcf, subplot(3,1, (mod((i-1),3) + 1)));
    hold on;
    %plot(T(1:length(TemplateMatchValue)),TemplateMatchValue/max(TemplateMatchValue) * 5000, Colours(1));
    %plot(T(1:length(TemplateMatchValue)),TemplateMatchValue2/max(TemplateMatchValue2) * 5000, Colours(7));
    %plot(T(1:length(TemplateMatchValue3)),TemplateMatchValue3/max(TemplateMatchValue3) * 5000, Colours(3));
    %plot(T(1:length(TemplateMatchValue)),Correlation/max(Correlation) * 5000, Colours(4));
    plot(T(1:length(TemplateMatchValue)),Correlation2/max(Correlation2) * 5000, Colours((mod((i-1), length(Colours)) + 1)));    
    %plot(T(1:length(TemplateMatchValue3)),Correlation3/max(Correlation3) * 5000, Colours(5));
    %plot(T(Indices(locs)),Correlation2(Indices(locs))/max(Correlation2) * 5000, [Colours(7),'s']);
end
