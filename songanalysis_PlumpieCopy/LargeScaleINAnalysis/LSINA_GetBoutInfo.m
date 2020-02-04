function [Bouts, Gaps, BoutDetails] = LSINA_GetBoutInfo(SURecordingDetails)

% This is done separately for triggered and continuous data. 
% I'm using a structure that keeps track of all the things related to a
% bout in different columns
% Column 1 and 2 are bout onset and offset indices
% Column 3 and 4 are the onset and offset files
% Column 5 and 6 are the onset and offset times of the bout with respect to
% time in that particular file, even if the data is continuous
% Column 7 is whether a motif is present (1) or not (0)
% Column 8 is >0 if there is >= 2000ms data in front, otherwise it is <0
% Column 9 is >1 if there is >= 2000ms data at the back, otherwise it is <1

% In addition, I am going to use the silent gaps for getting spontaneous
% activity - so I should keep track of gaps here - will do this later.
% For gaps:
% Column 1 - length of gap in ms
% Column 2 and 3 - onset and offset times of gap with respect to time in
% that particular file, even if the data is continuous
% Column 4 and 5 - onset and offset files nos
% Column 6 and 7 - onset and offset time for gap in terms of total time;
% this is only relevant for continuous data

% Finally, I also need to get rid of door open, close events when the
% female was put in or taken out. I also need to remove female calls from
% this

% 01.06.2017 - Added a new output called BoutDetails that will have the
% labels, onsets and offsets for each bout in a separate field. The onsets
% and offsets are relative to the first syllable onset of the bout

Bouts = [];
Gaps = [];
BoutDetails = [];

if (SURecordingDetails.Continuousdata == 0)
    for i = 1:length(SURecordingDetails.NoteInfo),
        if (isempty(SURecordingDetails.NoteInfo{i}.onsets))
            % If there are no onsets or offsets, then the entire file
            % constitutes a long silent gap and I should include this in
            % the gaps file.
            Gaps = [Gaps; [SURecordingDetails.FileLen(i) 0 SURecordingDetails.FileLen(i) i i]];
            continue;
        end
        % Calculate inter-syllable intervals
        Intervals = SURecordingDetails.NoteInfo{i}.onsets(2:end) - SURecordingDetails.NoteInfo{i}.offsets(1:end-1);
        
        % Now put in the data for gaps
        % While putting in the data, should remember that the silent period
        % before the first syllable is also a gap and the silent period
        % after the last syllable is also a gap
        
        if (~isempty(Intervals))
            TempOffsets = SURecordingDetails.NoteInfo{i}.offsets(1:end);
            TempOffsets = [0; TempOffsets(:)]; 
            TempOnsets = SURecordingDetails.NoteInfo{i}.onsets(1:end);
            TempOnsets = [TempOnsets(:); SURecordingDetails.FileLen(i)];
            
            Gaps = [Gaps; [[SURecordingDetails.NoteInfo{i}.onsets(1); Intervals(:); (SURecordingDetails.FileLen(i) - SURecordingDetails.NoteInfo{i}.offsets(end))] TempOffsets(:) TempOnsets(:) ones(size(TempOnsets(:)))*i ones(size(TempOffsets(:)))*i]];
        else
            % If intervals is empty, then it means that there was only 1
            % syllable in the entire file. Then, I can include all the time
            % before that syllable and all the time after that syllable as
            % a gap
            TempGapLengths = [SURecordingDetails.NoteInfo{i}.onsets(1); (SURecordingDetails.FileLen(i) - SURecordingDetails.NoteInfo{i}.offsets(1))];
            TempGapOnsets = [0; SURecordingDetails.NoteInfo{i}.offsets(1)];
            TempGapOffsets = [SURecordingDetails.NoteInfo{i}.onsets(1); SURecordingDetails.FileLen(i)];
            Gaps = [Gaps; [TempGapLengths(:) TempGapOnsets(:) TempGapOffsets(:) ones(size(TempGapOnsets))*i ones(size(TempGapOffsets))*i]];
        end
        
        % See if any are greater than the inter-bout interval
        LongIntervals = find(Intervals >= SURecordingDetails.Interboutinterval);
        
        if (isempty(LongIntervals))
            % This corresponds to all syllables in the file being part of
            % one long bout
            
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
            
            BoutDetails{end+1}.Labels = SURecordingDetails.NoteInfo{i}.labels(1:end);
            BoutDetails{end}.Onsets = SURecordingDetails.NoteInfo{i}.onsets(1:end) - SURecordingDetails.NoteInfo{i}.onsets(1);
            BoutDetails{end}.Offsets = SURecordingDetails.NoteInfo{i}.offsets(1:end) - SURecordingDetails.NoteInfo{i}.onsets(1);
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
            Bouts(end,8:9) = [(((Bouts(end,5) - SURecordingDetails.Interboutinterval)/SURecordingDetails.Interboutinterval)) (((SURecordingDetails.NoteInfo{i}.onsets(LongIntervals(1) + 1) - Bouts(end,6))/SURecordingDetails.Interboutinterval))];
            
            BoutDetails{end+1}.Labels = SURecordingDetails.NoteInfo{i}.labels(1:LongIntervals(1));
            BoutDetails{end}.Onsets = SURecordingDetails.NoteInfo{i}.onsets(1:LongIntervals(1)) - SURecordingDetails.NoteInfo{i}.onsets(1);
            BoutDetails{end}.Offsets = SURecordingDetails.NoteInfo{i}.offsets(1:LongIntervals(1)) - SURecordingDetails.NoteInfo{i}.onsets(1);
            
            for k = 2:length(LongIntervals),
                MotifFlag = 0; 
                for j = 1:length(SURecordingDetails.MotifLabels),
                    if (~isempty(find(SURecordingDetails.NoteInfo{i}.labels((LongIntervals(k-1) + 1):LongIntervals(k)) == SURecordingDetails.MotifLabels(j))))
                        MotifFlag = 1;
                        break;
                    end
                end

                Bouts(end+1,1:2) = [(LongIntervals(k-1) + 1) LongIntervals(k)];
                Bouts(end,3:4) = [i i];
                Bouts(end,5:6) = [SURecordingDetails.NoteInfo{i}.onsets(LongIntervals(k-1) + 1) SURecordingDetails.NoteInfo{i}.offsets(LongIntervals(k))];
                Bouts(end,7) = MotifFlag;
                Bouts(end,8:9) = [(((Bouts(end,5) - SURecordingDetails.Interboutinterval)/SURecordingDetails.Interboutinterval)) (((SURecordingDetails.NoteInfo{i}.onsets(LongIntervals(k) + 1) - Bouts(end,6))/SURecordingDetails.Interboutinterval))];
                
                BoutDetails{end+1}.Labels = SURecordingDetails.NoteInfo{i}.labels((LongIntervals(k-1) + 1):LongIntervals(k));
                BoutDetails{end}.Onsets = SURecordingDetails.NoteInfo{i}.onsets((LongIntervals(k-1) + 1):LongIntervals(k)) - SURecordingDetails.NoteInfo{i}.onsets((LongIntervals(k-1) + 1));
                BoutDetails{end}.Offsets = SURecordingDetails.NoteInfo{i}.offsets((LongIntervals(k-1) + 1):LongIntervals(k)) - SURecordingDetails.NoteInfo{i}.onsets((LongIntervals(k-1) + 1));
            end
            
            MotifFlag = 0; 
            for j = 1:length(SURecordingDetails.MotifLabels),
                if (~isempty(find(SURecordingDetails.NoteInfo{i}.labels((LongIntervals(end) + 1):end) == SURecordingDetails.MotifLabels(j))))
                    MotifFlag = 1;
                    break;
                end
            end
            
            Bouts(end+1,1:2) = [(LongIntervals(end) + 1) length(SURecordingDetails.NoteInfo{i}.labels)];
            Bouts(end,3:4) = [i i];
            Bouts(end,5:6) = [SURecordingDetails.NoteInfo{i}.onsets(LongIntervals(end) + 1) SURecordingDetails.NoteInfo{i}.offsets(end)];
            Bouts(end,7) = MotifFlag;
            Bouts(end,8:9) = [(((Bouts(end,5) - SURecordingDetails.Interboutinterval)/SURecordingDetails.Interboutinterval)) (((SURecordingDetails.FileLen(i) - Bouts(end,6))/SURecordingDetails.Interboutinterval))];
            
            BoutDetails{end+1}.Labels = SURecordingDetails.NoteInfo{i}.labels((LongIntervals(end) + 1):end);
            BoutDetails{end}.Onsets = SURecordingDetails.NoteInfo{i}.onsets((LongIntervals(end) + 1):end) - SURecordingDetails.NoteInfo{i}.onsets((LongIntervals(end) + 1));
            BoutDetails{end}.Offsets = SURecordingDetails.NoteInfo{i}.offsets((LongIntervals(end) + 1):end) - SURecordingDetails.NoteInfo{i}.onsets((LongIntervals(end) + 1));
        end
    end
else
    [AllLabels, AllOnsets, AllOffsets, AllUnAdjustedOnsets, AllUnAdjustedOffsets, OnsetFileNos, OffsetFileNos] = CombineContinuousDataNoteFiles(SURecordingDetails.DataDirectory, SURecordingDetails.SongFileNames, fullfile(SURecordingDetails.DataDirectory, 'ASSLNoteFiles'), SURecordingDetails.FileType);
    
    % Now get rid of female in, out events and female calls
    if (isfield(SURecordingDetails, 'Femalecalllabel'))
        FemaleINOUTEvents = regexp(AllLabels, ['[', SURecordingDetails.Femalecalllabel, SURecordingDetails.FemaleinLabel, SURecordingDetails.FemaleoutLabel, SURecordingDetails.Closedoorlabel, ']']);
        if (~isempty(FemaleINOUTEvents))
            AllLabels(FemaleINOUTEvents) = [];
            AllOnsets(FemaleINOUTEvents) = [];
            AllOffsets(FemaleINOUTEvents) = [];
            AllUnAdjustedOnsets(FemaleINOUTEvents) = [];
            AllUnAdjustedOffsets(FemaleINOUTEvents) = [];
            OnsetFileNos(FemaleINOUTEvents) = [];
            OffsetFileNos(FemaleINOUTEvents) = [];
        end
    end
    
    % Make everything into column vectors
    AllOnsets = AllOnsets(:);
    AllOffsets = AllOffsets(:);
    AllUnAdjustedOnsets = AllUnAdjustedOnsets(:);
    AllUnAdjustedOffsets = AllUnAdjustedOffsets(:);
    OnsetFileNos = OnsetFileNos(:);
    OffsetFileNos = OffsetFileNos(:);
    
    Intervals = AllOnsets(2:end) - AllOffsets(1:end-1);
    
    % Now fill in data for the gaps
    Gaps = [Gaps; [Intervals(:) AllUnAdjustedOffsets(1:end-1) AllUnAdjustedOnsets(2:end) OffsetFileNos(1:end-1) OnsetFileNos(2:end) AllOffsets(1:end-1) AllOnsets(2:end)]];
    
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
        Bouts(end,8:9) = [((AllOnsets(1) - SURecordingDetails.Interboutinterval)/SURecordingDetails.Interboutinterval) (((sum(SURecordingDetails.FileLen) - AllOffsets(end))/SURecordingDetails.Interboutinterval))];
        
        BoutDetails{end+1}.Labels = AllLabels(1:end);
        BoutDetails{end}.Onsets = AllOnsets(1:end) - AllOnsets(1);
        BoutDetails{end}.Offsets = AllOffsets(1:end) - AllOnsets(1);
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
        Bouts(end,8:9) = [(((AllOnsets(1) - SURecordingDetails.Interboutinterval)/SURecordingDetails.Interboutinterval)) ((AllOnsets(LongIntervals(1) + 1) - AllOffsets(LongIntervals(1)))/SURecordingDetails.Interboutinterval)];
        
        BoutDetails{end+1}.Labels = AllLabels(1:LongIntervals(1));
        BoutDetails{end}.Onsets = AllOnsets(1:LongIntervals(1)) - AllOnsets(1);
        BoutDetails{end}.Offsets = AllOffsets(1:LongIntervals(1)) - AllOnsets(1);
        
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
            Bouts(end,8:9) = [((Intervals(LongIntervals(k-1)) - SURecordingDetails.Interboutinterval)/SURecordingDetails.Interboutinterval) ((AllOnsets(LongIntervals(k) + 1) - AllOffsets(LongIntervals(k)))/SURecordingDetails.Interboutinterval)];
            
            BoutDetails{end+1}.Labels = AllLabels((LongIntervals(k-1) + 1):LongIntervals(k));
            BoutDetails{end}.Onsets = AllOnsets((LongIntervals(k-1) + 1):LongIntervals(k)) - AllOnsets((LongIntervals(k-1) + 1));
            BoutDetails{end}.Offsets = AllOffsets((LongIntervals(k-1) + 1):LongIntervals(k)) - AllOnsets((LongIntervals(k-1) + 1));
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
        Bouts(end,8:9) = [((Intervals(LongIntervals(end)) - SURecordingDetails.Interboutinterval)/SURecordingDetails.Interboutinterval) ((sum(SURecordingDetails.FileLen) - AllOffsets(end))/SURecordingDetails.Interboutinterval)];
        
        BoutDetails{end+1}.Labels = AllLabels((LongIntervals(end) + 1):end);
        BoutDetails{end}.Onsets = AllOnsets((LongIntervals(end) + 1):end) - AllOnsets((LongIntervals(end) + 1));
        BoutDetails{end}.Offsets = AllOffsets((LongIntervals(end) + 1):end) - AllOnsets((LongIntervals(end) + 1));
    end
end
