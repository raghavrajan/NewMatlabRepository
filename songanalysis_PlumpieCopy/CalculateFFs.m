function [] = CalculateFFs(DirectoryName, FileList)

Fid = fopen(FileList, 'r');
FileName = fgetl(Fid);
Index = 1;

cd(DirectoryName);

sm_win = 0.008;

FFTWinSize = sm_win; % in sec
FFTWinOverlap = 0.9;
Cols = ['rgbcmk'];

while (ischar(FileName(1)))
    load(FileName);
    if ((isfield(AnalysisOutput.DirFileInfo, 'MotifFiringRate')) && (isfield(AnalysisOutput.UnDirFileInfo, 'MotifFiringRate')))
        for j = 1:size(AnalysisOutput.UnDirFileInfo.WEventParameters.UWFirstSpikeTime,2),
            clear DirAmplitude UnDirAmplitude;
            for i = 1:size(AnalysisOutput.DirFileInfo.WEventParameters.UWFirstSpikeTime,1),
                if (~isempty(strfind(AnalysisOutput.FileType, 'okrank')))
                    [RawSong, Fs] = SSAReadOKrankData(AnalysisOutput.RawDataDirectory, AnalysisOutput.RecFileDirectory, AnalysisOutput.DirFileInfo.FileNames{AnalysisOutput.DirFileInfo.Syllables.Index(i)},1);
                else
                    if (~isempty(strfind(AnalysisOutput.FileType, 'obs')))
                        [RawSong, Fs] = SSASoundIn(AnalysisOutput.RawDataDirectory, AnalysisOutput.RecFileDirectory, AnalysisOutput.DirFileInfo.FileNames{AnalysisOutput.DirFileInfo.Syllables.Index(i)}, 'obs0r');
                        RawSong = RawSong * 5/32768;
                    end
                end
                Time = (1:1:length(RawSong))/Fs;
                %BurstOnsetTime = AnalysisOutput.DirFileInfo.WEventParameters.UWFirstSpikeTime(i,j) + AnalysisOutput.DirFileInfo.Syllables.Start(i,1);
                BurstOnsetTime = mean(AnalysisOutput.DirFileInfo.WEventParameters.FirstSpikeTime(:,j)) + AnalysisOutput.DirFileInfo.Syllables.Start(i,1);
                BurstSong = RawSong(find((Time >= (BurstOnsetTime - 0.1)) & (Time <= (BurstOnsetTime + 0.1))));

                WinSize = round(FFTWinSize * Fs);
                WinOverlap = round(FFTWinOverlap * WinSize);
                [S, F, T, P] = spectrogram(BurstSong, hamming(WinSize), WinOverlap, WinSize, Fs);

                Freq1 = find((F >= 860) & (F <= 8600));
                DirAmplitude(i,:) = log10(sum(S(Freq1,:).*conj(S(Freq1,:))));
            end
            
            for i = 1:size(AnalysisOutput.UnDirFileInfo.WEventParameters.UWFirstSpikeTime,1),
                if (~isempty(strfind(AnalysisOutput.FileType, 'okrank')))
                    [RawSong, Fs] = SSAReadOKrankData(AnalysisOutput.RawDataDirectory, AnalysisOutput.RecFileDirectory, AnalysisOutput.UnDirFileInfo.FileNames{AnalysisOutput.UnDirFileInfo.Syllables.Index(i)},1);
                else
                    if (~isempty(strfind(AnalysisOutput.FileType, 'obs')))
                        [RawSong, Fs] = SSASoundIn(AnalysisOutput.RawDataDirectory, AnalysisOutput.RecFileDirectory, AnalysisOutput.UnDirFileInfo.FileNames{AnalysisOutput.UnDirFileInfo.Syllables.Index(i)}, 'obs0r');
                        RawSong = RawSong * 5/32768;
                    end
                end
                Time = (1:1:length(RawSong))/Fs;
                %BurstOnsetTime = AnalysisOutput.UnDirFileInfo.WEventParameters.UWFirstSpikeTime(i,j) + AnalysisOutput.UnDirFileInfo.Syllables.Start(i,1);
                BurstOnsetTime = mean(AnalysisOutput.UnDirFileInfo.WEventParameters.FirstSpikeTime(:,j)) + AnalysisOutput.UnDirFileInfo.Syllables.Start(i,1);
                BurstSong = RawSong(find((Time >= (BurstOnsetTime - 0.1)) & (Time <= (BurstOnsetTime + 0.1))));

                WinSize = round(FFTWinSize * Fs);
                WinOverlap = round(FFTWinOverlap * WinSize);
                [S, F, T, P] = spectrogram(BurstSong, hamming(WinSize), WinOverlap, WinSize, Fs);

                Freq1 = find((F >= 860) & (F <= 8600));
                UnDirAmplitude(i,:) = log10(sum(S(Freq1,:).*conj(S(Freq1,:))));
            end

            DirT = linspace(-0.1, 0.1, size(DirAmplitude, 2));
            UnDirT = linspace(-0.1, 0.1, size(UnDirAmplitude, 2));
            MinSpikes = min([min(AnalysisOutput.DirFileInfo.WEventParameters.NoofSpikes(:,j)) min(AnalysisOutput.UnDirFileInfo.WEventParameters.NoofSpikes(:,j))]);
            MaxSpikes = max([max(AnalysisOutput.DirFileInfo.WEventParameters.NoofSpikes(:,j)) max(AnalysisOutput.UnDirFileInfo.WEventParameters.NoofSpikes(:,j))]);

            FigNo = gcf;
            for i = MinSpikes:1:MaxSpikes,
                figure(FigNo + 1);
                subplot(2,2,1);
                hold on;
                if (length(find(AnalysisOutput.DirFileInfo.WEventParameters.NoofSpikes(:,j) == i)) > 0)
                    plot(DirT, DirAmplitude(find(AnalysisOutput.DirFileInfo.WEventParameters.NoofSpikes(:,j) == i), :)', Cols(i-MinSpikes+1));
                end
            end
            subplot(2,2,2);
            for i = MinSpikes:1:MaxSpikes,
                figure(FigNo + 1);
                hold on;
                if (length(find(AnalysisOutput.DirFileInfo.WEventParameters.NoofSpikes(:,j) == i)) > 0)
                   plot(DirT, mean(DirAmplitude(find(AnalysisOutput.DirFileInfo.WEventParameters.NoofSpikes(:,j) == i), :))', Cols(i-MinSpikes+1));
                end
            end
            subplot(2,2,3);
            for i = MinSpikes:1:MaxSpikes,
                figure(FigNo + 1);
                hold on;
                if (length(find(AnalysisOutput.UnDirFileInfo.WEventParameters.NoofSpikes(:,j) == i)) > 0)
                    plot(UnDirT, UnDirAmplitude(find(AnalysisOutput.UnDirFileInfo.WEventParameters.NoofSpikes(:,j) == i), :)', Cols(i-MinSpikes+1));
                end
            end
            subplot(2,2,4);
            for i = MinSpikes:1:MaxSpikes,
                figure(FigNo + 1);
                hold on;
                if (length(find(AnalysisOutput.UnDirFileInfo.WEventParameters.NoofSpikes(:,j) == i)) > 0)
                    plot(UnDirT, mean(UnDirAmplitude(find(AnalysisOutput.UnDirFileInfo.WEventParameters.NoofSpikes(:,j) == i), :))', Cols(i-MinSpikes+1));
                end
            end
            subplot(2,2,1);
            title([FileName, ' Burst #', num2str(j), ' aligned on warped mean first spike time']);
        end
    end
    FileName = fgetl(Fid);
end

fclose(Fid);