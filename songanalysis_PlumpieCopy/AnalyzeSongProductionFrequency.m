function [] = AnalyzeSongProductionFrequency(DirectoryName, FileList, FileType)

PresentDir = pwd;

Fid = fopen(FileList, 'r');
Temp = textscan(Fid, '%s', 'delimiter', '\n');
fclose(Fid);

SongFiles = Temp{1};

FFTWinSize = 8; % in seconds
FFTWinOverlap = 0.5;

AllAmplitudes = [];
for SongFileNo = 1:length(SongFiles),
    SongFile = SongFiles{SongFileNo};

    if (isempty(SongFile))
        continue;
    end
    
    Slash = find((SongFile == '/') | (SongFile == '\'));
    if (~isempty(Slash))
        SongFile = SongFile(Slash(end)+1:end);
    end
    disp(SongFile);

    try
        if (strfind(FileType, 'okrank'))
            [Song, Fs] = ReadOKrankData(DirectoryName, SongFile, 1);
        else
            if (strfind(FileType, 'wav'))
                cd(DirectoryName);
                [Song, Fs] = wavread(SongFile);
                cd(PresentDir);
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
    [LogAmplitude] = ASSLCalculateLogAmplitude(Song, Fs, [1:1:length(Song)]/Fs, FFTWinSize, FFTWinOverlap);
    LogAmplitude = resample(LogAmplitude, 100, Fs);
    AllAmplitudes = [AllAmplitudes; LogAmplitude(:)];
end

SongFreq = fft(AllAmplitudes);
figure;
plot(100/2*linspace(0, 1, length(AllAmplitudes)/2+1), 2*abs(SongFreq(1:length(AllAmplitudes)/2+1)));

disp('Finished');


    