function [] = PlotMakeTemplatesMatch(DirectoryName, BatchFile, FFTBinSize, FileType, ChanNo, Template)

if (DirectoryName(end) ~= '/')
    DirectoryName(end + 1) = '/';
end

cd(DirectoryName);

NoteInfo = [];
fid = fopen(BatchFile,'r');

FileNames = textscan(fid, '%s', 'DeLimiter', '\n');
FileNames = FileNames{1};

fclose(fid);

for i = 1:length(FileNames),
    FileName = FileNames{i};
    BirdName = FileName(1:(find((FileName == '_') | (FileName == '-'))));
    if (strfind(FileType, 'okrank'))
        [RawData, Fs] = ReadOKrankData(DirectoryName, FileName, ChanNo);
    else
        if (strfind(FileType, 'obs'))
            [RawData, Fs] = soundin_copy(DirectoryName, FileName, 'obs0r');    
        end
    end
    Time = (1:1:length(RawData))/Fs;
    TempIndices = 1:(FFTBinSize*Fs/1000):length(RawData)-(FFTBinSize*Fs/1000);
    clear Match MatchTime;
    for k = 1:length(TempIndices);
        TempData = RawData((round((k-1)*FFTBinSize*Fs/1000) + 1):(round((k)*FFTBinSize*Fs/1000)));
        TempFFT = abs(fft(TempData));
        TempFFT(1) = 0;
        %TempFFT = TempFFT(1:(length(TempFFT)/2+1))/sum(TempFFT(1:(length(TempFFT)/2+1)));
        TempFFT = TempFFT(1:(length(TempFFT)/2+1));
        TempFFT = (TempFFT - mean(TempFFT))/std(TempFFT);
        Match(k) = length(Template)/sum((TempFFT - Template).*(TempFFT - Template));
        MatchTime(k) = Time(TempIndices(k));
    end
    figure;
    plot(MatchTime, Match);
    PlotSpectrogram(DirectoryName, FileName, FileType);
    hold on;
    plot(MatchTime, ((Match - mean(Match))/std(Match))*500 + 4000, 'b');
    uiwait(gcf);
end

disp('Finished');