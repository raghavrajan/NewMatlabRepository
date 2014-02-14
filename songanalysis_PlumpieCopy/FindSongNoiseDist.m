function [Threshold] = FindSongNoiseDist(DirectoryName, FileList, FileType)

PresentDirectory = pwd;

Fid = fopen(FileList, 'r');
SongFile = fgetl(Fid);

cd(DirectoryName);

Edges = -5:0.001:5;
SongFileIndex = 1;

while (ischar(SongFile(1)))
    Slash = find((SongFile == '\') | (SongFile == '/'));
    
    if (~isempty(Slash))
        SongFile = SongFile(Slash(end)+1:end);
    end
    disp(SongFile);
    try
        if (strfind(FileType, 'okrank'))
            [Song, Fs] = ReadOKrankData(DirectoryName, SongFile, 1);
        else
            if (strfind(FileType, 'wav'))
                [Song, Fs] = wavread(SongFile);
            else
                if (strfind(FileType, 'obs'))
                    channel_string = strcat('obs',num2str(0),'r');
                    [Song, Fs] = soundin_copy(DirectoryName, SongFile, channel_string);
    
                    % Convert to uV - 5V on the data acquisition is 32768
                    Song = Song * 5/32768;
                end
            end
        end
    catch
        continue;
    end
    if (isempty(Song))
        continue;
    end
    SongTime = (1:1:length(Song))/Fs;
    
    FFTWinSize = 0.008; % in sec
    FFTWinOverlap = 0.9;
    
    WinSize = round(FFTWinSize * Fs);
    WinOverlap = round(FFTWinOverlap * WinSize);

    [S, F, T, P] = spectrogram(Song, hamming(WinSize), WinOverlap, WinSize, Fs);
    
    Freq1 = find((F >= 860) & (F <= 8600));
    Power = log10(sum(S(Freq1,:).*conj(S(Freq1,:))));
    Power = spline(T, Power, SongTime);
    Obj = gmdistribution.fit(Power', 2);
    [Sound(SongFileIndex,1), Index] = max(Obj.mu);
    Sound(SongFileIndex,2) = Obj.Sigma(Index);
    [Noise(SongFileIndex,1), Index] = min(Obj.mu);
    Noise(SongFileIndex,2) = Obj.Sigma(Index);
    SongFileIndex = SongFileIndex + 1;
    SongFile = fgetl(Fid);
end

fclose(Fid);
cd(PresentDirectory);