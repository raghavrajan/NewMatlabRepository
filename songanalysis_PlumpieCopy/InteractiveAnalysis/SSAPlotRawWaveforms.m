function [RawWaveformFigure] = SSAPlotRawWaveforms(DirectoryName, RecFileDirectoryName, DirFileInfo, UnDirFileInfo, ChannelNo, FileType, PreSongStartDuration, NoofExamples)

RawWaveformFigure = figure;
axes('Position',[0.15 0.15 0.75 0.75]);
set(gca,'Box','off');
hold on;
MaxRawData = 0;
set(gcf, 'Color', 'w');

if (~isempty(DirFileInfo))
    [DirFileInfo.SongsLengths, SongIndices] = sort(DirFileInfo.SongLengths);
    DirFileInfo.Syllables.Start = DirFileInfo.Syllables.Start(SongIndices, :);
    DirFileInfo.Syllables.End = DirFileInfo.Syllables.End(SongIndices, :);
    DirFileInfo.Syllables.Index = DirFileInfo.Syllables.Index(SongIndices);
end

if (~isempty(DirFileInfo))
    if (size(DirFileInfo.Syllables.Start,1) > NoofExamples)
        NoofSongs = NoofExamples;
    else
        NoofSongs = size(DirFileInfo.Syllables.Start,1);
    end
        
    for i = 1:NoofSongs,
        FileIndex = DirFileInfo.Syllables.Index(i);

        if (strfind(FileType, 'obs'))
            [RawData,Fs] = SSASoundIn(DirectoryName, RecFileDirectoryName, [DirFileInfo.FileNames{FileIndex}],['obs',num2str(DirFileInfo.NChans{FileIndex} - ChannelNo),'r']);
            RawData = RawData * 500/32768;
        else
            if (strfind(FileType, 'okrank'))
                [RawData,Fs] = SSAReadOKrankData(DirectoryName, RecFileDirectoryName, [DirFileInfo.FileNames{FileIndex}], num2str(ChannelNo));
                RawData = RawData * 100;
            end
        end

        [b, a] = butter(4, [300/16000 6000/16000]);
        RawData = filtfilt(b, a, RawData);
        Time = (1:1:length(RawData))/Fs;
        Indices = find((Time > (DirFileInfo.Syllables.Start(i,1) - PreSongStartDuration)) & (Time < (DirFileInfo.Syllables.End(i,end))));
        TempMax = max(RawData(Indices)) - min(RawData(Indices));
        TempMax = TempMax + 0.05 * TempMax;
        RawData = RawData + MaxRawData;
        MaxRawData = MaxRawData + TempMax;

        plot((Time(Indices) - Time(Indices(1)) - PreSongStartDuration), RawData(Indices),'r');

        if (isfield(DirFileInfo, 'SpikeData'));
            SpikeIndices = find((DirFileInfo.SpikeData.Times{FileIndex} > (DirFileInfo.Syllables.Start(i,1) - PreSongStartDuration)) & (DirFileInfo.SpikeData.Times{FileIndex} < (DirFileInfo.Syllables.End(i,end))));
        
            for j = 1:length(SpikeIndices),
                SpikeTimeIndices = find((Time <= (DirFileInfo.SpikeData.Times{FileIndex}(SpikeIndices(j)))), 1, 'last');
                plot((Time((SpikeTimeIndices-8):(SpikeTimeIndices+23)) - Time(Indices(1)) - PreSongStartDuration), RawData((SpikeTimeIndices-8):(SpikeTimeIndices+23)), 'b'); 
            end
        end
    end
end

if (~isempty(UnDirFileInfo))
    [UnDirFileInfo.SongsLengths, SongIndices] = sort(UnDirFileInfo.SongLengths);
    UnDirFileInfo.Syllables.Start = UnDirFileInfo.Syllables.Start(SongIndices, :);
    UnDirFileInfo.Syllables.End = UnDirFileInfo.Syllables.End(SongIndices, :);
    UnDirFileInfo.Syllables.Index = UnDirFileInfo.Syllables.Index(SongIndices);
end

if (~isempty(UnDirFileInfo))
    if (size(UnDirFileInfo.Syllables.Start,1) > NoofExamples)
        NoofSongs = NoofExamples;
    else
        NoofSongs = size(UnDirFileInfo.Syllables.Start,1);
    end
        
    for i = 1:NoofSongs,
        FileIndex = UnDirFileInfo.Syllables.Index(i);

        if (strfind(FileType, 'obs'))
            [RawData,Fs] = SSASoundIn(DirectoryName, RecFileDirectoryName, [UnDirFileInfo.FileNames{FileIndex}],['obs',num2str(UnDirFileInfo.NChans{FileIndex} - ChannelNo),'r']);
            RawData = RawData * 500/32768;
        else
            if (strfind(FileType, 'okrank'))
                [RawData,Fs] = SSAReadOKrankData(DirectoryName, RecFileDirectoryName, [UnDirFileInfo.FileNames{FileIndex}],num2str(ChannelNo));
                RawData = RawData * 100;
            end
        end
        [b, a] = butter(4, [300/16000 6000/16000]);
        RawData = filtfilt(b, a, RawData);
        
        Time = (1:1:length(RawData))/Fs;
        Indices = find((Time > (UnDirFileInfo.Syllables.Start(i,1) - PreSongStartDuration)) & (Time < (UnDirFileInfo.Syllables.End(i,end))));
        TempMax = max(RawData(Indices)) - min(RawData(Indices));
        TempMax = TempMax + 0.05 * TempMax;
        RawData = RawData + MaxRawData;
        MaxRawData = MaxRawData + TempMax;

        plot((Time(Indices) - Time(Indices(1)) - PreSongStartDuration), RawData(Indices),'k');
        
        if (isfield(UnDirFileInfo, 'SpikeData'));
            SpikeIndices = find((UnDirFileInfo.SpikeData.Times{FileIndex} > (UnDirFileInfo.Syllables.Start(i,1) - PreSongStartDuration)) & (UnDirFileInfo.SpikeData.Times{FileIndex} < (UnDirFileInfo.Syllables.End(i,end))));

            for j = 1:length(SpikeIndices),
                SpikeTimeIndices = find((Time <= (UnDirFileInfo.SpikeData.Times{FileIndex}(SpikeIndices(j)))), 1, 'last');
                plot((Time((SpikeTimeIndices-8):(SpikeTimeIndices+23)) - Time(Indices(1)) - PreSongStartDuration), RawData((SpikeTimeIndices-8):(SpikeTimeIndices+23)), 'b'); 
            end
        end
    end
end

if (~isempty(DirFileInfo))
    if (strfind(FileType, 'observer'))
        DotIndex = find(DirFileInfo.FileNames{1} == '.');
        TitleName = DirFileInfo.FileNames{1}(1:(DotIndex(1) - 1));
    else
        if (strfind(FileType, 'okrank'))
            TitleName = DirFileInfo.FileNames{1}(1:(end - 6));
        end
    end
else
    if (~isempty(UnDirFileInfo))
        if (strfind(FileType, 'observer'))
            DotIndex = find(UnDirFileInfo.FileNames{1} == '.');
            TitleName = UnDirFileInfo.FileNames{1}(1:(DotIndex(1) - 1));
        else
            if (strfind(FileType, 'okrank'))
                TitleName = UnDirFileInfo.FileNames{1}(1:(end - 6));
            end
        end
    end
end

axis tight;
set(gca,'FontSize',14,'FontWeight','bold');
xlabel('Time (sec)','FontSize',14,'FontWeight','bold');
ylabel('Amplitude (uV)','FontSize',14,'FontWeight','bold');
title([TitleName,' - Raw Waveforms'],'FontSize',16,'FontWeight','bold');        
    