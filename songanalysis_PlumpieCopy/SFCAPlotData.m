 function [handles] = SFCAPlotData(handles)

% First plot the song
cla(handles.SFCA_Axes);
PlotSpectrogramInAxisLowRes(handles.SFCA.DirName, handles.SFCA.SongFiles{handles.SFCA.SongFilesIndex}, handles.SFCA.FileTypes{handles.SFCA.ChosenFileType}, handles.SFCA_Axes);
xticks('auto');
xlabel('Time (sec)');
yticks([]);
ylabel([]);
handles.SFCA.PlotAxisLimits = axis;

% Now add the notes
cla(handles.SFCA_LabelAxes);
axes(handles.SFCA_LabelAxes);
for i = 1:length(handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.onsets),
    text((handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.onsets(i))/1000, 0, handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.labels(i), 'FontSize', 14);
end
handles.SFCA.LabelAxisLimits = handles.SFCA.PlotAxisLimits;
handles.SFCA.LabelAxisLimits(3:4) = [-0.5 0.5];
axis(handles.SFCA.LabelAxisLimits);
xticks([]);
yticks([]);
title(handles.SFCA.SongFiles{handles.SFCA.SongFilesIndex}, 'FontSize', 14);

% Now to check if motif labels are present and make them red
if (~isempty(handles.SFCA.MotifLabels))
    Motifs = strfind(handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.labels, handles.SFCA.MotifLabels);
    if (~isempty(Motifs))
        MotifOnsets = handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.onsets(Motifs)/1000;
        MotifOffsets = handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.offsets(Motifs + length(handles.SFCA.MotifLabels) - 1)/1000;
        
        for i = 1:length(Motifs),
            axes(handles.SFCA_Axes);
            hold on;
            patch([MotifOnsets(i)*ones(1,2) MotifOffsets(i)*ones(1,2)], [handles.SFCA.PlotAxisLimits(3:4) fliplr(handles.SFCA.PlotAxisLimits(3:4))], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
            
            axes(handles.SFCA_LabelAxes);
            hold on;
            patch([MotifOnsets(i)*ones(1,2) MotifOffsets(i)*ones(1,2)], [handles.SFCA.LabelAxisLimits(3:4) fliplr(handles.SFCA.LabelAxisLimits(3:4))], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        end
    end
end

Bouts = [];
% Now check for any intervals greater than the inter-bout interval and draw
% a patch around these in blue
Intervals = (handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.onsets(2:end) - handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.offsets(1:end-1))/1000;
LongIntervals = find(Intervals >= handles.SFCA.InterBoutInterval);
if (~isempty(LongIntervals))
    for i = 1:length(LongIntervals),
        if (i == 1)
            if (handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.onsets(1)/1000 >= handles.SFCA.InterBoutInterval)
                Bouts(end+1,:) = [1 LongIntervals(i)];
            end
        else
            Bouts(end+1,:) = [LongIntervals(i-1)+1 LongIntervals(i)];
        end
        axes(handles.SFCA_Axes);
        hold on;
        PatchXVals = [handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.offsets(LongIntervals(i))/1000*ones(1,2) handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.onsets(LongIntervals(i)+1)/1000*ones(1,2)];
        patch(PatchXVals, [handles.SFCA.PlotAxisLimits(3:4) fliplr(handles.SFCA.PlotAxisLimits(3:4))], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.1);

        axes(handles.SFCA_LabelAxes);
        hold on;
        patch(PatchXVals, [handles.SFCA.LabelAxisLimits(3:4) fliplr(handles.SFCA.LabelAxisLimits(3:4))], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
    end
    Bouts(end+1,:) = [LongIntervals(end)+1 length(handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.labels)];
else
    Bouts(end+1,:) = [1 length(handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.labels)];
end

if (handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.onsets(1)/1000 >= handles.SFCA.InterBoutInterval)
    axes(handles.SFCA_Axes);
    hold on;
    PatchXVals = [0*ones(1,2) handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.onsets(1)/1000*ones(1,2)];
    patch(PatchXVals, [handles.SFCA.PlotAxisLimits(3:4) fliplr(handles.SFCA.PlotAxisLimits(3:4))], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
    
    axes(handles.SFCA_LabelAxes);
    hold on;
    patch(PatchXVals, [handles.SFCA.LabelAxisLimits(3:4) fliplr(handles.SFCA.LabelAxisLimits(3:4))], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
end

if (handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.offsets(end)/1000 <= (handles.SFCA.FileLen(handles.SFCA.SongFilesIndex) - handles.SFCA.InterBoutInterval))
    axes(handles.SFCA_Axes);
    hold on;
    PatchXVals = [handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.offsets(end)/1000*ones(1,2) handles.SFCA.FileLen(handles.SFCA.SongFilesIndex)*ones(1,2)];
    patch(PatchXVals, [handles.SFCA.PlotAxisLimits(3:4) fliplr(handles.SFCA.PlotAxisLimits(3:4))], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
    
    axes(handles.SFCA_LabelAxes);
    hold on;
    patch(PatchXVals, [handles.SFCA.LabelAxisLimits(3:4) fliplr(handles.SFCA.LabelAxisLimits(3:4))], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
end

% Now to check if IN labels are present and make them green - I want to
% make the entire IN period green based on the last few INs that are <
% 500ms apart - have to do this for each new bout

for i = 1:size(Bouts,1),
    FirstMotifSyll = NaN;
    if ((~isempty(handles.SFCA.INLabels)) && (~isempty(handles.SFCA.MotifLabels)))
        % First find first motif syllable
        for j = Bouts(i,1):Bouts(i,2),
            if (~isempty(find(handles.SFCA.MotifLabels == handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.labels(j))))
                FirstMotifSyll = j;
                break;
            end
        end

        % Now find INs
        if isempty(find(isnan(FirstMotifSyll)))
            INs = zeros(1,FirstMotifSyll-Bouts(i,1));
            for j = Bouts(i,1):FirstMotifSyll-1,
                if (~isempty(find(handles.SFCA.INLabels == handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.labels(j))))
                    INs(j) = 1;
                end
            end
            NotINs = find(INs == 0);
            if (~isempty(NotINs))
                NotINs = NotINs(end);
                if (NotINs == length(INs))
                    LastINs = [];
                else
                    LastINs = NotINs+1:1:length(INs);
                end
            else
                LastINs = 1:1:length(INs);
            end
            % Now find the last set of INs with < 500ms between them
            if (length(LastINs) >= 2)
                IntervalsBetweenLastINs = handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.onsets(LastINs(2:end)) - handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.offsets(LastINs(1:end-1));
                LongINIntervals = find(IntervalsBetweenLastINs >= handles.SFCA.InterINInterval*1000);
                if (~isempty(LongINIntervals))
                    LastINs = LastINs(LongINIntervals+1:end);
                end
            end
            % Now to draw the green patch
            if (~isempty(LastINs))
                LastINOnsets = handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.onsets(LastINs(1))/1000;
                LastINOffsets = handles.SFCA.NoteInfo{handles.SFCA.SongFilesIndex}.offsets(LastINs(end))/1000;

                axes(handles.SFCA_Axes);
                hold on;
                patch([LastINOnsets*ones(1,2) LastINOffsets*ones(1,2)], [handles.SFCA.PlotAxisLimits(3:4) fliplr(handles.SFCA.PlotAxisLimits(3:4))], 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.3);

                axes(handles.SFCA_LabelAxes);
                hold on;
                patch([LastINOnsets*ones(1,2) LastINOffsets*ones(1,2)], [handles.SFCA.LabelAxisLimits(3:4) fliplr(handles.SFCA.LabelAxisLimits(3:4))], 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
            end
        end
    end
end
