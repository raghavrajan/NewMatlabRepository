function [SeqMatchValues] = CalculateMotifSeqMatch(DirectoryName, FileList, FileType, MotifTemplate, PlotOption)

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
    Freq = find((F>=860) & (F<=8600));
    %Data = ScaleSpect(S);
    S = log10(abs(S(Freq,:)));
    disp(S(1));
    [Match] = CalTemplateMatch(MotifTemplate, S);
%     for i = 1:(size(S,2) - size(MotifTemplate, 2) + 1),
%         WinIndices = i:1:(i + size(MotifTemplate, 2) - 1);
%         Data = ScaleSpect(S);
%         %Data = log10(abs(S(Freq, WinIndices)));
%          Data = S(:, WinIndices);
%         Data = (Data - mean(mean(Data)))/std(reshape(Data, 1, (size(Data,1) * size(Data,2))));
%         Match(i,1) = 1/sum(sum(abs(Data - MotifTemplate)));
%        if (i == 1)
%             disp(['Mean is ', num2str(mean(mean(Data))), ' and std is ', num2str(std(reshape(Data, 1, (size(Data,1)*size(Data,2)))))]);
%         end
%     end
    SeqMatchValues(end + 1) = max(Match);
    disp(['SongFile is ', SongFile, ' and max sequence match value is ', num2str(SeqMatchValues(end))]);
    if (strfind(PlotOption, 'on'))
        PlotSpectrogram(DirectoryName, SongFile, FileType);
        hold on;
        plot(T(1:length(Match)), (Match/max(Match) *8000), 'r');
        zoom xon;
        uiwait(gcf);
    end
    clear Match;
    SongFile = fgetl(Fid);
end
fclose(Fid);