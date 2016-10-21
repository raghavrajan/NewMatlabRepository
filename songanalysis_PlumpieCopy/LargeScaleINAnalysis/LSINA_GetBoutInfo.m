function [Bouts, Gaps] = LSINA_GetBoutInfo(SURecordingDetails)

% This is done separately for triggered and continuous data. 
% I'm using a structure that keeps track of all the things related to a
% bout in different columns
% Column 1 and 2 are bout onset and offset indices
% Column 3 and 4 are the onset and offset files
% Column 5 and 6 are the onset and offset times of the bout with respect to
% time in that particular file, even if the data is continuous
% Column 7 is whether a motif is present (1) or not (0)
% Column 8 is >1 if there is >= 2000ms d7:9ata in front, otherwise it is 0
% Column 9 is >1 if there is >= 2000ms data at the back, otherwise it is 0

% In addition, I am going to use the silent gaps for getting spontaneous
% activity - so I should keep track of gaps here - will do this later.

Bouts = [];

if (SURecordingDetails.Continuousdata == 0)
    for i = 1:length(SURecordingDetails.NoteInfo),
        if (isempty(SURecordingDetails.NoteInfo{i}.onsets))
            continue;
        end
        % Calculate inter-syllable intervals
        Intervals = SURecordingDetails.NoteInfo{i}.onsets(2:end) - SURecordingDetails.NoteInfo{i}.offsets(1:end-1);
        
        % See if any are greater than the inter-bout interval
        LongIntervals = find(Intervals >= SURecordingDetails.Interboutinterval);
        
        if (isempty(LongIntervals))
            MotifFlag = 0; 
            for j = 1:length(SURecordingDetails.MotifLabels),
                if (~isempty(find(SURecordingDetails.NoteInfo{i}.labels == SURecordingDetails.MotifLabels(j))))
                    MotifFlag = 1;
                    break;
                end
            end
            
            Bouts(end+1,1:2) = [1 length(SURecordingDetails.NoteInfo{i}.labels)];
            Bouts(end,3:4) = [i i];
            Bouts(end,5:6) = [SURecordingDetails.NoteInfo{i}.onsets(1) SURecordingDetails.NoteInfo{i}.offsets(end)];
            Bouts(end,7) = MotifFlag;
            Bouts(end,8:9) = [(((Bouts(end,5) - SURecordingDetails.Interboutinterval)/SURecordingDetails.Interboutinterval)) (((SURecordingDetails.FileLen(i) - Bouts(end,6))/SURecordingDetails.Interboutinterval))];
        else
            MotifFlag = 0; 
            for j = 1:length(SURecordingDetails.MotifLabels),
                if (~isempty(find(SURecordingDetails.NoteInfo{i}.labels(1:LongIntervals(1)) == SURecordingDetails.MotifLabels(j))))
                    MotifFlag = 1;
                    break;
                end
            end
            
            Bouts(end+1,1:2) = [1 LongIntervals(1)];
            Bouts(end,3:4) = [i i];
            Bouts(end,5:6) = [SURecordingDetails.NoteInfo{i}.onsets(1) SURecordingDetails.NoteInfo{i}.offsets(LongIntervals(1))];
            Bouts(end,7) = MotifFlag;
            Bouts(end,8:9) = [(((Bouts(end,5) - SURecordingDetails.Interboutinterval)/SURecordingDetails.Interboutinterval)) (((SURecordingDetails.FileLen(i) - Bouts(end,6))/SURecordingDetails.Interboutinterval))];
            
            for k = 2:length(LongIntervals),
                MotifFlag = 0; 
                for j = 1:length(SURecordingDetails.MotifLabels),
                    if (~isempty(find(SURecordingDetails.NoteInfo{i}.labels(LongIntervals(k-1)+1:LongIntervals(k)) == SURecordingDetails.MotifLabels(j))))
                        MotifFlag = 1;
                        break;
                    end
                end

                Bouts(end+1,1:2) = [(LongIntervals(k-1) + 1) LongIntervals(k)];
                Bouts(end,3:4) = [i i];
                Bouts(end,5:6) = [SURecordingDetails.NoteInfo{i}.onsets(LongIntervals(k-1) + 1) SURecordingDetails.NoteInfo{i}.offsets(LongIntervals(k))];
                Bouts(end,7) = MotifFlag;
                Bouts(end,8:9) = [(((Bouts(end,5) - SURecordingDetails.Interboutinterval)/SURecordingDetails.Interboutinterval)) (((SURecordingDetails.FileLen(i) - Bouts(end,6))/SURecordingDetails.Interboutinterval))];
            end
            
            MotifFlag = 0; 
            for j = 1:length(SURecordingDetails.MotifLabels),
                if (~isempty(find(SURecordingDetails.NoteInfo{i}.labels(LongIntervals(end)+1:end) == SURecordingDetails.MotifLabels(j))))
                    MotifFlag = 1;
                    break;
                end
            end
            
            Bouts(end+1,1:2) = [(LongIntervals(end) + 1) length(SURecordingDetails.NoteInfo{i}.labels)];
            Bouts(end,3:4) = [i i];
            Bouts(end,5:6) = [SURecordingDetails.NoteInfo{i}.onsets(LongIntervals(end) + 1) SURecordingDetails.NoteInfo{i}.offsets(end)];
            Bouts(end,7) = MotifFlag;
            Bouts(end,8:9) = [(((Bouts(end,5) - SURecordingDetails.Interboutinterval)/SURecordingDetails.Interboutinterval)) (((SURecordingDetails.FileLen(i) - Bouts(end,6))/SURecordingDetails.Interboutinterval))];
        end
    end
else
    [AllLabels, AllOnsets, AllOffsets, AllUnAdjustedOnsets, AllUnAdjustedOffsets, OnsetFileNos, OffsetFileNos] = CombineContinuousDataNoteFiles(SURecordingDetails.DataDirectory, SURecordingDetails.SongFileNames, fullfile(SURecordingDetails.DataDirectory, 'ASSLNoteFiles'), SURecordingDetails.FileType);
    
    Intervals = AllOnsets(2:end) - AllOffsets(1:end-1);
        
    % See if any are greater than the inter-bout interval
    LongIntervals = find(Intervals >= SURecordingDetails.Interboutinterval);
        
    if (isempty(LongIntervals))
        MotifFlag = 0; 
        for j = 1:length(SURecordingDetails.MotifLabels),
            if (~isempty(find(AllLabels == SURecordingDetails.MotifLabels(j))))
                MotifFlag = 1;
                break;
            end
        end

        Bouts(end+1,1:2) = [1 length(AllLabels)];
        Bouts(end,3:4) = [OnsetFileNos(1) OffsetFileNos(end)];
        Bouts(end,5:6) = [AllUnAdjustedOnsets(1) AllUnAdjustedOffsets(end)];
        Bouts(end,7) = MotifFlag;
        Bouts(end,8:9) = [(((Bouts(end,5) - SURecordingDetails.Interboutinterval)/SURecordingDetails.Interboutinterval)) (((SURecordingDetails.FileLen(OffsetFileNos(end)) - Bouts(end,6))/SURecordingDetails.Interboutinterval))];
    else
        MotifFlag = 0; 
        for j = 1:length(SURecordingDetails.MotifLabels),
            if (~isempty(find(AllLabels(1:LongIntervals(1)) == SURecordingDetails.MotifLabels(j))))
                MotifFlag = 1;
                break;
            end
        end

        Bouts(end+1,1:2) = [1 LongIntervals(1)];
        Bouts(end,3:4) = [OnsetFileNos(1) OffsetFileNos(LongIntervals(1))];
        Bouts(end,5:6) = [AllUnAdjustedOnsets(1) AllUnAdjustedOffsets(LongIntervals(1))];
        Bouts(end,7) = MotifFlag;
        if (OnsetFileNos(1) == 1)
            Bouts(end,8:9) = [(((Bouts(end,5) - SURecordingDetails.Interboutinterval)/SURecordingDetails.Interboutinterval)) 1];
        else
            Bouts(end,8:9) = [1 1];
        end
        
        for k = 2:length(LongIntervals),
            MotifFlag = 0; 
            for j = 1:length(SURecordingDetails.MotifLabels),
                if (~isempty(find(AllLabels(LongIntervals(k-1)+1:LongIntervals(k)) == SURecordingDetails.MotifLabels(j))))
                    MotifFlag = 1;
                    break;
                end
            end

            Bouts(end+1,1:2) = [(LongIntervals(k-1) + 1) LongIntervals(k)];
            Bouts(end,3:4) = [OnsetFileNos(LongIntervals(k-1) + 1) OffsetFileNos(LongIntervals(k))];
            Bouts(end,5:6) = [AllUnAdjustedOnsets(LongIntervals(k-1) + 1) AllUnAdjustedOffsets(LongIntervals(k))];
            Bouts(end,7) = MotifFlag;
            Bouts(end,8:9) = [1 1];
        end

        MotifFlag = 0; 
        for j = 1:length(SURecordingDetails.MotifLabels),
            if (~isempty(find(AllLabels(LongIntervals(end)+1:end) == SURecordingDetails.MotifLabels(j))))
                MotifFlag = 1;
                break;
            end
        end

        Bouts(end+1,1:2) = [(LongIntervals(end) + 1) length(AllLabels)];
        Bouts(end,3:4) = [OnsetFileNos(LongIntervals(end) + 1) OffsetFileNos(end)];
        Bouts(end,5:6) = [AllUnAdjustedOnsets(LongIntervals(end) + 1) AllUnAdjustedOffsets(end)];
        Bouts(end,7) = MotifFlag;
        
        if (OffsetFileNos(end) == length(SURecordingDetails.SongFileNames))
            Bouts(end,8:9) = [1 (((SURecordingDetails.FileLen(OffsetFileNos(end)) - Bouts(end,6))/SURecordingDetails.Interboutinterval))];
        else
            Bouts(end,8:9) = [1 1];
        end
    end
   
end
