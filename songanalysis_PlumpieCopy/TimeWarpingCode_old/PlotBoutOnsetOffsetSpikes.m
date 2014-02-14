function [] = PlotBoutOnsetOffsetSpikes(DirFileInfo,UnDirFileInfo,Motif,MainFigure);

% First do the bout onsets and offsets for directed bouts

AllNoteOnsets = [0; [DirFileInfo.Notes.NoteOnsets]];
BoutBeginIndices = find(diff(AllNoteOnsets) > 4);
if (length(BoutBeginIndices) > 0)
    BoutBeginTimes = DirFileInfo.Notes.NoteOnsets(BoutBeginIndices);
    DirBoutBeginSpikeTimes = [];
    for i = 1:length(BoutBeginTimes)
        TempSpikeTimes = DirFileInfo.SpikeData.Times(find(([DirFileInfo.SpikeData.Times] > (BoutBeginTimes(i) - 4)) & ([DirFileInfo.SpikeData.Times] < BoutBeginTimes(i)))) - BoutBeginTimes(i);
        DirBoutBeginSpikeTimes = [DirBoutBeginSpikeTimes; [TempSpikeTimes ones(length(TempSpikeTimes),1)*i]];
    end
else
    DirBoutBeginSpikeTimes = [];
end

AllNoteOffsets = [[DirFileInfo.Notes.NoteOffsets]; sum(DirFileInfo.RecordLengths)];
BoutEndIndices = find(diff(AllNoteOffsets) > 4);
if (length(BoutEndIndices) > 0)
    BoutEndTimes = DirFileInfo.Notes.NoteOffsets(BoutEndIndices);
    DirBoutEndSpikeTimes = [];
    for i = 1:length(BoutEndTimes)
        TempSpikeTimes = DirFileInfo.SpikeData.Times(find(([DirFileInfo.SpikeData.Times] > (BoutEndTimes(i))) & ([DirFileInfo.SpikeData.Times] < (BoutEndTimes(i) + 4)))) - BoutEndTimes(i);
        DirBoutEndSpikeTimes = [DirBoutEndSpikeTimes; [TempSpikeTimes ones(length(TempSpikeTimes),1)*i]];
    end
else
    DirBoutEndSpikeTimes = [];    
end

AllNoteOnsets = [0; [UnDirFileInfo.Notes.NoteOnsets]];
BoutBeginIndices = find(diff(AllNoteOnsets) > 4);
if (length(BoutBeginIndices) > 0)
    BoutBeginTimes = UnDirFileInfo.Notes.NoteOnsets(BoutBeginIndices);
    UnDirBoutBeginSpikeTimes = [];
    for i = 1:length(BoutBeginTimes)
        TempSpikeTimes = UnDirFileInfo.SpikeData.Times(find(([UnDirFileInfo.SpikeData.Times] > (BoutBeginTimes(i) - 4)) & ([UnDirFileInfo.SpikeData.Times] < BoutBeginTimes(i)))) - BoutBeginTimes(i);
        UnDirBoutBeginSpikeTimes = [UnDirBoutBeginSpikeTimes; [TempSpikeTimes ones(length(TempSpikeTimes),1)*i]];
    end
else
    UnDirBoutBeginSpikeTimes = [];
end

if (length(DirBoutBeginSpikeTimes) > 0)
    UnDirBoutBeginSpikeTimes(:,2) = UnDirBoutBeginSpikeTimes(:,2) + DirBoutBeginSpikeTimes(end,2) + 2;
end

AllNoteOffsets = [[UnDirFileInfo.Notes.NoteOffsets]; sum(UnDirFileInfo.RecordLengths)];
BoutEndIndices = find(diff(AllNoteOffsets) > 4);
if (length(BoutEndIndices) > 0)
    BoutEndTimes = UnDirFileInfo.Notes.NoteOffsets(BoutEndIndices);
    UnDirBoutEndSpikeTimes = [];
    for i = 1:length(BoutEndTimes)
        TempSpikeTimes = UnDirFileInfo.SpikeData.Times(find(([UnDirFileInfo.SpikeData.Times] > (BoutEndTimes(i))) & ([UnDirFileInfo.SpikeData.Times] < (BoutEndTimes(i) + 4)))) - BoutEndTimes(i);
        UnDirBoutEndSpikeTimes = [UnDirBoutEndSpikeTimes; [TempSpikeTimes ones(length(TempSpikeTimes),1)*i]];
    end
else
    UnDirBoutEndSpikeTimes = [];
end

if (length(DirBoutEndSpikeTimes) > 0)
    UnDirBoutEndSpikeTimes(:,2) = UnDirBoutEndSpikeTimes(:,2) + DirBoutEndSpikeTimes(end,2) + 2;
end



figure(MainFigure);
BoutOnsetPlot = axes('Position',[0.15 0.15 0.75 0.35]);
set(gca,'Box','off');
hold on;

BoutOffsetPlot = axes('Position',[0.15 0.55 0.75 0.35]);
set(gca,'Box','off');
hold on;

axes(BoutOnsetPlot);

if (length(DirBoutBeginSpikeTimes) > 0)
    plot(DirBoutBeginSpikeTimes(:,1),DirBoutBeginSpikeTimes(:,2),'w+');
    MarkerString = repmat('|',size(DirBoutBeginSpikeTimes(:,1),1),1);
    text(DirBoutBeginSpikeTimes(:,1),DirBoutBeginSpikeTimes(:,2),MarkerString,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',2,'FontName','FixedWidth','FontUnits','pixels','Color','k','FontWeight','bold');
end

if (length(UnDirBoutBeginSpikeTimes) > 0)
    plot(UnDirBoutBeginSpikeTimes(:,1),UnDirBoutBeginSpikeTimes(:,2),'w+');
    MarkerString = repmat('|',size(UnDirBoutBeginSpikeTimes(:,1),1),1);
    text(UnDirBoutBeginSpikeTimes(:,1),UnDirBoutBeginSpikeTimes(:,2),MarkerString,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',2,'FontName','FixedWidth','FontUnits','pixels','Color',[0.6 0.6 0.6],'FontWeight','bold');
end

axis([-4 0 0 UnDirBoutBeginSpikeTimes(end,2)]);

set(gca,'FontSize',14,'FontWeight','bold');

xlabel('Time (sec)','FontSize',16,'FontWeight','bold');
ylabel('Trials','FontSize',16,'FontWeight','bold');

axes(BoutOffsetPlot);

if (length(DirBoutEndSpikeTimes) > 0)
    plot(DirBoutEndSpikeTimes(:,1),DirBoutEndSpikeTimes(:,2),'w+');
    MarkerString = repmat('|',size(DirBoutEndSpikeTimes(:,1),1),1);
    text(DirBoutEndSpikeTimes(:,1),DirBoutEndSpikeTimes(:,2),MarkerString,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',2,'FontName','FixedWidth','FontUnits','pixels','Color','k','FontWeight','bold');
end

if (length(UnDirBoutEndSpikeTimes) > 0)
    plot(UnDirBoutEndSpikeTimes(:,1),UnDirBoutEndSpikeTimes(:,2),'w+');
    MarkerString = repmat('|',size(UnDirBoutEndSpikeTimes(:,1),1),1);
    text(UnDirBoutEndSpikeTimes(:,1),UnDirBoutEndSpikeTimes(:,2),MarkerString,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',2,'FontName','FixedWidth','FontUnits','pixels','Color',[0.6 0.6 0.6],'FontWeight','bold');
end

axis([0 4 0 UnDirBoutEndSpikeTimes(end,2)]);

set(gca,'FontSize',14,'FontWeight','bold');

ylabel('Trials','FontSize',16,'FontWeight','bold');

disp('Finished doing bout onsets and offsets');