function [] = PlotEventWaveforms(DirectoryName, ChannelNo, DirFileInfo, UnDirFileInfo, MedianMotif, MainFigure)

figure(MainFigure);

if (isfield(DirFileInfo.EventParameters,'EventTimes') > 0)
    EventPlotWidth = (0.8 - 0.02*(length(DirFileInfo.EventParameters.EventTimes) - 1))/length(DirFileInfo.EventParameters.EventTimes);
    NoofEvents = length(DirFileInfo.EventParameters.EventTimes);
else
    if (isfield(UnDirFileInfo.EventParameters,'EventTimes') > 0)
        NoofEvents = length(UnDirFileInfo.EventParameters.EventTimes);
        EventPlotWidth = (0.8 - 0.02*(length(DirFileInfo.EventParameters.EventTimes) - 1))/length(UnDirFileInfo.EventParameters.EventTimes);
    else
        return;
    end
end

for EventCounter = 1:NoofEvents,
    axes('Position',[((EventCounter-1)*(EventPlotWidth + 0.02) + 0.15) 0.4 EventPlotWidth 0.35]);
    set(gca,'Box','off');
    hold on;
    MaxRawData = 0;

    if (length(DirFileInfo.RecordLengths) > 0)
        TempEventTime = DirFileInfo.EventParameters.EventTimes(EventCounter);
        if (size(DirFileInfo.Syllables.Start,1) > 2)
            NoofSongs = 2;
        else
            NoofSongs = size(DirFileInfo.Syllables.Start,1);
        end
        
        for i = 1:NoofSongs,
            FileIndex = find(cumsum(DirFileInfo.RecordLengths) < DirFileInfo.Syllables.Start(i,1),1,'last');
    
            if (length(FileIndex) == 0)
                FileIndex = 1;
                TotalRecordLength = 0;
           else
                FileIndex = FileIndex + 1;
                TotalRecordLength = sum(DirFileInfo.RecordLengths(1:(FileIndex - 1)));
            end

            [RawData,Fs] = soundin_copy(DirectoryName, [DirFileInfo.FileNames{FileIndex}],['obs',num2str(ChannelNo),'r']);

            RawData = RawData * 500/32768;

            Time = 0:1/Fs:length(RawData)/Fs;
            
            MotifPosition = find(MedianMotif.SyllableStartings > TempEventTime,1,'first');
            if (TempEventTime > (MedianMotif.GapStartings(MotifPosition - 1)))
                WarpedPosition = (TempEventTime - MedianMotif.GapStartings(MotifPosition - 1))/(MedianMotif.SyllableStartings(MotifPosition) - MedianMotif.GapStartings(MotifPosition - 1));
                WarpedPosition = WarpedPosition * DirFileInfo.Gaps.Length(i,(MotifPosition - 1));
                StartTime = DirFileInfo.Gaps.Start(i,(MotifPosition - 1)) + WarpedPosition - 0.015;
                EndTime = DirFileInfo.Gaps.Start(i,(MotifPosition - 1)) + WarpedPosition + 0.015;
            else
                WarpedPosition = (TempEventTime - MedianMotif.SyllableStartings(MotifPosition - 1))/(MedianMotif.GapStartings(MotifPosition - 1) - MedianMotif.SyllableStartings(MotifPosition - 1));
                WarpedPosition = WarpedPosition * DirFileInfo.Syllables.Length(i,(MotifPosition - 1));
                StartTime = DirFileInfo.Syllables.Start(i,(MotifPosition - 1)) + WarpedPosition - 0.015 - TotalRecordLength;
                EndTime = DirFileInfo.Syllables.Start(i,(MotifPosition - 1)) + WarpedPosition + 0.015 - TotalRecordLength;
            end
            
            Indices = find((Time > StartTime) & (Time < EndTime));                
            
            TempMax = max(RawData(Indices)) - min(RawData(Indices));
            TempMax = TempMax + 0.05 * TempMax;
            RawData = RawData + MaxRawData;
            MaxRawData = MaxRawData + TempMax;

            plot((Time(Indices) - Time(Indices(1))),RawData(Indices),'k');
        end
    end

    if (length(UnDirFileInfo.RecordLengths) > 0)
        TempEventTime = UnDirFileInfo.EventParameters.EventTimes(EventCounter);
       for i = 1:2,
            FileIndex = find(cumsum(UnDirFileInfo.RecordLengths) < UnDirFileInfo.Syllables.Start(i,1),1,'last');

            if (length(FileIndex) == 0)
                FileIndex = 1;
                TotalRecordLength = 0;
            else
                FileIndex = FileIndex + 1;
                TotalRecordLength = sum(UnDirFileInfo.RecordLengths(1:(FileIndex - 1)));
            end

            [RawData,Fs] = soundin_copy(DirectoryName, [UnDirFileInfo.FileNames{FileIndex}],['obs',num2str(ChannelNo),'r']);

            RawData = RawData * 500/32768;

            Time = 0:1/Fs:length(RawData)/Fs;
            MotifPosition = find(MedianMotif.SyllableStartings > TempEventTime,1,'first');
            if (TempEventTime > (MedianMotif.GapStartings(MotifPosition - 1)))
                WarpedPosition = (TempEventTime - MedianMotif.GapStartings(MotifPosition - 1))/(MedianMotif.SyllableStartings(MotifPosition) - MedianMotif.GapStartings(MotifPosition - 1));
                WarpedPosition = WarpedPosition * DirFileInfo.Gaps.Length(i,(MotifPosition - 1));
                StartTime = UnDirFileInfo.Gaps.Start(i,(MotifPosition - 1)) + WarpedPosition - 0.015;
                EndTime = UnDirFileInfo.Gaps.Start(i,(MotifPosition - 1)) + WarpedPosition + 0.015;
            else
                WarpedPosition = (TempEventTime - MedianMotif.SyllableStartings(MotifPosition - 1))/(MedianMotif.GapStartings(MotifPosition - 1) - MedianMotif.SyllableStartings(MotifPosition - 1));
                WarpedPosition = WarpedPosition * DirFileInfo.Syllables.Length(i,(MotifPosition - 1));
                StartTime = UnDirFileInfo.Syllables.Start(i,(MotifPosition - 1)) + WarpedPosition - 0.015 - TotalRecordLength;
                EndTime = UnDirFileInfo.Syllables.Start(i,(MotifPosition - 1)) + WarpedPosition + 0.015 - TotalRecordLength;
            end
            
            Indices = find((Time > StartTime) & (Time < EndTime));                

            TempMax = max(RawData(Indices)) - min(RawData(Indices));
            TempMax = TempMax + 0.05 * TempMax;
            RawData = RawData + MaxRawData;
            MaxRawData = MaxRawData + TempMax;

            plot((Time(Indices) - Time(Indices(1))), RawData(Indices),'Color',[0.6 0.6 0.6]);
            hold on;
            axis tight; 
            set(gca,'Visible','off');
        end
    end
end

