function [] = PlotBoutSpikeTimes(DirFileInfo,UnDirFileInfo,Motif,MainFigure)

figure(MainFigure);

RasterPlot = axes('Position',[0.1 0.1 0.85 0.85]);
hold on;

DirFileLengths = cumsum(DirFileInfo.RecordLengths);
UnDirFileLengths = cumsum(UnDirFileInfo.RecordLengths);

RasterYValue = 0;

if (length(DirFileInfo.PST) > 0)
    for i = 1:length(DirFileInfo.FileNames),
        if (i == 1)
            TempSpikeTrain = DirFileInfo.SpikeData.Times(find((DirFileInfo.SpikeData.Times > 0) & (DirFileInfo.SpikeData.Times <= DirFileLengths(i))));
            TempNoteLabels = DirFileInfo.Notes.NoteLabels(find((DirFileInfo.Notes.NoteOnsets > 0) & (DirFileInfo.Notes.NoteOnsets <= DirFileLengths(i))));
            TempNoteOnsets = DirFileInfo.Notes.NoteOnsets(find((DirFileInfo.Notes.NoteOnsets > 0) & (DirFileInfo.Notes.NoteOnsets <= DirFileLengths(i))));
        else
            TempSpikeTrain = DirFileInfo.SpikeData.Times(find((DirFileInfo.SpikeData.Times > DirFileLengths(i-1)) & (DirFileInfo.SpikeData.Times <= DirFileLengths(i))));
            TempNoteLabels = DirFileInfo.Notes.NoteLabels(find((DirFileInfo.Notes.NoteOnsets > DirFileLengths(i-1)) & (DirFileInfo.Notes.NoteOnsets <= DirFileLengths(i))));            
            TempNoteOnsets = DirFileInfo.Notes.NoteOnsets(find((DirFileInfo.Notes.NoteOnsets > DirFileLengths(i-1)) & (DirFileInfo.Notes.NoteOnsets <= DirFileLengths(i))));            
            if (length(TempNoteOnsets) > 0)
                TempNoteOnsets = TempNoteOnsets - DirFileLengths(i-1);
            end
            
            if (length(TempSpikeTrain) > 0)
                TempSpikeTrain = TempSpikeTrain - DirFileLengths(i-1);
            end
            
        end

        if (length(TempNoteOnsets) > 0)
            if (strfind(TempNoteLabels,Motif))
                for i = 1:length(TempNoteOnsets),
                    plot(TempNoteOnsets(i), 1*RasterYValue,'w+');
                    text(TempNoteOnsets(i),1*RasterYValue,TempNoteLabels(i),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',2,'FontName','FixedWidth','FontUnits','pixels','Color','k','FontWeight','bold');
                end
                RasterYValue = RasterYValue + 1;
                if (length(TempSpikeTrain) > 0)
                    plot(TempSpikeTrain, ones(length(TempSpikeTrain),1)*RasterYValue,'w+');
                    MarkerString = repmat('|',length(TempSpikeTrain),1);
                    text(TempSpikeTrain,ones(length(TempSpikeTrain),1)*RasterYValue,MarkerString,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',2,'FontName','FixedWidth','FontUnits','pixels','Color','k','FontWeight','bold');
                    RasterYValue = RasterYValue + 1;
                end
            end
        end
    end
end

if (length(UnDirFileInfo.PST) > 0)
    for i = 1:length(UnDirFileInfo.FileNames),
        if (i == 1)
            TempSpikeTrain = UnDirFileInfo.SpikeData.Times(find((UnDirFileInfo.SpikeData.Times > 0) & (UnDirFileInfo.SpikeData.Times <= UnDirFileLengths(i))));
            TempNoteLabels = UnDirFileInfo.Notes.NoteLabels(find((UnDirFileInfo.Notes.NoteOnsets > 0) & (UnDirFileInfo.Notes.NoteOnsets <= UnDirFileLengths(i))));
            TempNoteOnsets = UnDirFileInfo.Notes.NoteOnsets(find((UnDirFileInfo.Notes.NoteOnsets > 0) & (UnDirFileInfo.Notes.NoteOnsets <= UnDirFileLengths(i))));
        else
            TempSpikeTrain = UnDirFileInfo.SpikeData.Times(find((UnDirFileInfo.SpikeData.Times > UnDirFileLengths(i-1)) & (UnDirFileInfo.SpikeData.Times <= UnDirFileLengths(i))));
            TempNoteLabels = UnDirFileInfo.Notes.NoteLabels(find((UnDirFileInfo.Notes.NoteOnsets > UnDirFileLengths(i-1)) & (UnDirFileInfo.Notes.NoteOnsets <= UnDirFileLengths(i))));            
            TempNoteOnsets = UnDirFileInfo.Notes.NoteOnsets(find((UnDirFileInfo.Notes.NoteOnsets > UnDirFileLengths(i-1)) & (UnDirFileInfo.Notes.NoteOnsets <= UnDirFileLengths(i))));            
            if (length(TempNoteOnsets) > 0)
                TempNoteOnsets = TempNoteOnsets - UnDirFileLengths(i-1);
            end
            
            if (length(TempSpikeTrain) > 0)
                TempSpikeTrain = TempSpikeTrain - UnDirFileLengths(i-1);
            end
            
        end

        if (length(TempNoteOnsets) > 0)
            if (strfind(TempNoteLabels,Motif))
                for i = 1:length(TempNoteOnsets),
                    plot(TempNoteOnsets(i), 1*RasterYValue,'w+');
                    text(TempNoteOnsets(i),1*RasterYValue,TempNoteLabels(i),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',2,'FontName','FixedWidth','FontUnits','pixels','Color','k','FontWeight','bold');
                end
                RasterYValue = RasterYValue + 1;
                if (length(TempSpikeTrain) > 0)
                    plot(TempSpikeTrain, ones(length(TempSpikeTrain),1)*RasterYValue,'w+');
                    MarkerString = repmat('|',length(TempSpikeTrain),1);
                    text(TempSpikeTrain,ones(length(TempSpikeTrain),1)*RasterYValue,MarkerString,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',2,'FontName','FixedWidth','FontUnits','pixels','Color','k','FontWeight','bold');
                    RasterYValue = RasterYValue + 1;
                end
            end
        end
    end
end