function [SyllableData] = GetSyllableListInfo(BirdParameters)

%========== Help text for GetSyllableListInfo.m ===========================
% Description: This is a function used to get a comprehensive list of
% syllables from a filelist, with the corresponding file #, onsets and
% offsets for each of those syllables. 
% Usage: [SyllableListInfo] = GetSyllableListInfo(BirdParameters)
% Inputs: 
%   1) BirdParameters is a structure that is used by a different program.
%   The important fields of this structure that are used by this script are
%   NoteInfo, FileLen, Continuousdata. NoteFileInfo is a cell array, one for each
%   file and has the contents of the corresponding note file. FileLen is an
%   array with the length of each file in seconds. Continuousdata is 1 for
%   continuous data and 0 for triggered data.
% Outputs:
%   1) SyllableData - has one row for each syllable with the first
%   column having syllable labels (as a number - char of this will give the
%   actual syllable label), second column has file # with onset of syllable 
%   , third column has file # with offset of syllable, fourth column has
%   syllable onset in ms with respect to start of file with onset and fifth
%   column has syllable offset in ms with respect to start of file with
%   offset
%   If the data is continuous, then this SyllableData matrix has two more
%   columns, for the onsets and offsets relative to the entire time of the
%   continuous recording.

% =========================================================================

SyllableData = [];
if (BirdParameters.Continuousdata == 0) 
    for i = 1:length(BirdParameters.NoteInfo),
        TempSyllableData = [double(BirdParameters.NoteInfo{i}.labels(:)) ones(length(BirdParameters.NoteInfo{i}.labels),1)*i ones(length(BirdParameters.NoteInfo{i}.labels),1)*i BirdParameters.NoteInfo{i}.onsets(:) BirdParameters.NoteInfo{i}.offsets(:)];
        SyllableData = [SyllableData; TempSyllableData];
    end
else
    ContinuousFileTime = 0;
    ContinuousFirstSyllFlag = 0;
    ContinuousLastSyllFlag = 0;
    for i = 1:length(BirdParameters.NoteInfo),
        % Now check if the last syllable is a capital letter. If it is then
        % check if there is a next file and if there is a next file, check
        % if the first syllable in the next file is also a capital letter -
        % if this is the case, then these two have to be merged. If a
        % syllable has to be merged, in the Syllable Data matrix, it will
        % have one file for its onset file in column 2 and the next file
        % for its offset file in column 3.
        if (~isempty(BirdParameters.NoteInfo{i}.labels))
            if (regexp(BirdParameters.NoteInfo{i}.labels(end), '[A-Z]'))
                if ((i < length(BirdParameters.NoteInfo)) && (~isempty(BirdParameters.NoteInfo{i+1}.labels)))
                    if (regexp(BirdParameters.NoteInfo{i+1}.labels(1), '[A-Z]'))
                        if (BirdParameters.NoteInfo{i}.labels(end) == BirdParameters.NoteInfo{i+1}.labels(1))
                            ContinuousLastSyllFlag = 1;
                        else
                            ContinuousLastSyllFlag = 0;
                        end
                    else
                        ContinuousLastSyllFlag = 0;
                    end
                else
                    ContinuousLastSyllFlag = 0;
                end
            else
                ContinuousLastSyllFlag = 0;
            end
        
            if ((ContinuousLastSyllFlag == 0) && (ContinuousFirstSyllFlag == 0))
                TempLabels = BirdParameters.NoteInfo{i}.labels';
                TempOnsets = BirdParameters.NoteInfo{i}.onsets; 
                TempOffsets = BirdParameters.NoteInfo{i}.offsets;
                TempOnsetFileNo = ones(length(TempLabels), 1)*i;
                TempOffsetFileNo = ones(length(TempLabels),1)*i;            
            else
                if ((ContinuousLastSyllFlag == 0))
                    TempLabels = BirdParameters.NoteInfo{i}.labels(2:end)';
                    TempOnsets = BirdParameters.NoteInfo{i}.onsets(2:end); 
                    TempOffsets = BirdParameters.NoteInfo{i}.offsets(2:end);
                    TempOnsetFileNo = ones(length(TempLabels), 1)*i;
                    TempOffsetFileNo = ones(length(TempLabels), 1)*i;            
                else
                    if (ContinuousFirstSyllFlag == 0)
                        TempLabels = [BirdParameters.NoteInfo{i}.labels(1:end-1)'; lower(BirdParameters.NoteInfo{i}.labels(end))];
                        TempOnsets = BirdParameters.NoteInfo{i}.onsets; 
                        TempOffsets = [BirdParameters.NoteInfo{i}.offsets(1:end-1); BirdParameters.NoteInfo{i+1}.offsets(1)];
                        TempOnsetFileNo = ones(length(TempLabels), 1)*i;
                        TempOffsetFileNo = [ones(length(TempLabels)-1, 1)*i; i+1];            
                    else                  
                        TempLabels = [BirdParameters.NoteInfo{i}.labels(2:end-1)'; lower(BirdParameters.NoteInfo{i}.labels(end))] ;
                        TempOnsets = BirdParameters.NoteInfo{i}.onsets(2:end); 
                        TempOffsets = [BirdParameters.NoteInfo{i}.offsets(2:end-1); BirdParameters.NoteInfo{i+1}.offsets(1)];
                        TempOnsetFileNo = [ones(length(TempLabels), 1)*i];
                        TempOffsetFileNo = [ones(length(TempLabels)-1, 1)*i; i+1];            
                    end
                end
            end
            TempFullDurOnsets = TempOnsets + ContinuousFileTime;
            if (ContinuousLastSyllFlag == 1)
                TempFullDurOffsets = TempOffsets + ContinuousFileTime;
                TempFullDurOffsets(end) = TempFullDurOffsets(end) + BirdParameters.FileLen(i);
            else
                TempFullDurOffsets = TempOffsets + ContinuousFileTime;
            end
            ContinuousFileTime = ContinuousFileTime + BirdParameters.FileLen(i);
            
            SyllableData = [SyllableData; [double(TempLabels) TempOnsetFileNo(:) TempOffsetFileNo(:) TempOnsets(:) TempOffsets(:) TempFullDurOnsets(:) TempFullDurOffsets(:)]];
            clear TempLabels TempOnsetFileNo TempOffsetFileNo TempOnsets TempOffsets
            if (ContinuousLastSyllFlag == 1)
                ContinuousFirstSyllFlag = 1;
            else
                ContinuousFirstSyllFlag = 0;
            end
        end
    end
end