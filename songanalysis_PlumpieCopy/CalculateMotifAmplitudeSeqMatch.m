function [SeqMatchValues] = CalculateMotifAmplitudeSeqMatch(DirectoryName, FileList, FileType, MotifTemplate, PlotOption)

%cd(DirectoryName);
Fid = fopen(FileList, 'r');
SongFile = fgetl(Fid);

SeqMatchValues = [];
while(ischar(SongFile(1)))
    if (ispc)
        Slash = find(SongFile == '\');
    else
        Slash = find(SongFile == '/');
    end
    if (~isempty(Slash))
        SongFile = SongFile(Slash(end)+1:end);
    end
    disp(SongFile);
    % Get Song Data from observer file
    if (strfind(FileType,'obs'))
        channel_string = strcat('obs',num2str(0),'r');
        [rawsong,Fs] = soundin_copy(DirectoryName, SongFile,channel_string);

        % Convert to uV - 5V on the data acquisition is 32768
        rawsong = rawsong * 5/32768;
    else
        if (strfind(FileType,'wav'));
            [rawsong, Fs] = wavread(SongFile);
        else 
            if (strfind(FileType, 'okrank'))
                [rawsong, Fs] = ReadOKrankData(DirectoryName, SongFile, 1);
            end
        end
    end
    
    FFTWinSize = 0.008;
    FFTWinOverlap = 0.9;
    WinSize = round(FFTWinSize * Fs);
    WinOverlap = round(FFTWinOverlap * WinSize);
    [S, F, T, P] = spectrogram(rawsong, hamming(WinSize), WinOverlap, WinSize, Fs);
    Freq = find((F>=2000) & (F<=7000));
    Power = sum(log10(abs(S(Freq,:))));
    
    WinIndices = zeros(length(MotifTemplate), (length(Power) - length(MotifTemplate) + 1));
    WinIndices(1,:) = 1:1:(length(Power) - length(MotifTemplate) + 1);
    for i = 2:length(MotifTemplate),
        WinIndices(i,:) = WinIndices(i-1,:) + 1;
    end
    Data = Power(WinIndices);
    SumData = repmat(sum(Data), size(Data, 1), 1);
    STDData = repmat(std(Data), size(Data, 1), 1);
    Data = (Data - SumData)./STDData;
    Template = repmat(MotifTemplate, 1, size(Data,2));
    Match = 1./sum(abs(Data - Template));
    SeqMatchValues(end + 1) = max(Match);
    if (strfind(PlotOption, 'on'))
        PlotSpectrogram(DirectoryName, SongFile, FileType);
        hold on;
        plot(T(1:length(Match)), (Match/max(Match) *8000), 'r');
        zoom xon;
        uiwait(gcf);
    end
    disp(['SongFile is ', SongFile, ' and max sequence match value is ', num2str(SeqMatchValues(end))]);
    SongFile = fgetl(Fid);
end
fclose(Fid);