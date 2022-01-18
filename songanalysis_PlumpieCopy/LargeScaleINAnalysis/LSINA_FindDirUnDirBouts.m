function [UnDirBoutDetails, DirBoutDetails] = LSINA_FindDirUnDirBouts(SURecordingDetails)

InterBoutInterval = SURecordingDetails.Interboutinterval;

% If this is continuous data, then put all spikes together in one long
% array.
if (SURecordingDetails.Continuousdata == 1)
    % Put all syllable labels together in one long array
    [AllLabels, AllOnsets, AllOffsets] = CombineContinuousDataNoteFiles(SURecordingDetails.DataDirectory, SURecordingDetails.SongFileNames, fullfile(SURecordingDetails.DataDirectory, 'ASSLNoteFiles'), SURecordingDetails.FileType);
     % Now get rid of female in, out events and female calls
    FemaleINOUTEvents = regexp(AllLabels, ['[', SURecordingDetails.Femalecalllabel, SURecordingDetails.FemaleinLabel, SURecordingDetails.FemaleoutLabel, SURecordingDetails.Closedoorlabel, ']']);
    if (~isempty(FemaleINOUTEvents))
        AllLabels(FemaleINOUTEvents) = [];
        AllOnsets(FemaleINOUTEvents) = [];
        AllOffsets(FemaleINOUTEvents) = [];
    end
end

% First find all undir bouts that have more than the inter-bout
% interval before the start and at the end
% 25/02/2017 - realised a few days ago that the 9th column of Bouts has to
% be > 1 for the bout to considered as having enough post-data - have to
% redo some of the analysis and figures because of this.

UnDirBouts = find((SURecordingDetails.BoutDirUnDir == 0) & (SURecordingDetails.Bouts(:,8) > 0)' & (SURecordingDetails.Bouts(:,9) > 1)');

if (isempty(UnDirBouts))
   UnDirBoutDetails = [];
end

% Now for each of the neurons, take the bouts with enough data at the
% beginning and find all the spikes before bout initiation in two different
% windows - -2 to -1.5s (spont. activity window) and -0.6 to -0.1s
% (pre-song window)
% In addition get the spike train, the raster and the pst for that bout.
% Keep track of bout # also

if (SURecordingDetails.Continuousdata == 1)
    Index = 1;

    for j = UnDirBouts(:)',
        
        UnDirBoutDetails.BoutIndex(Index) = j;
        
        UnDirBoutDetails.BoutLabels{Index} = AllLabels(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        UnDirBoutDetails.BoutOnsets{Index} = AllOnsets(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        UnDirBoutDetails.BoutOffsets{Index} = AllOffsets(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        UnDirBoutDetails.BoutFirstSyllLabel(Index) = UnDirBoutDetails.BoutLabels{Index}(1);
        UnDirBoutDetails.SongBoutOrNot(Index) = SURecordingDetails.Bouts(j,7);
        
        Index = Index + 1;
    end
else
    Index = 1;
    for j = UnDirBouts(:)',
        UnDirBoutDetails.BoutIndex(Index) = j;
        
        UnDirBoutDetails.BoutLabels{Index} = SURecordingDetails.NoteInfo{SURecordingDetails.Bouts(j,3)}.labels(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        UnDirBoutDetails.BoutOnsets{Index} = SURecordingDetails.NoteInfo{SURecordingDetails.Bouts(j,3)}.onsets(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        UnDirBoutDetails.BoutOffsets{Index} = SURecordingDetails.NoteInfo{SURecordingDetails.Bouts(j,3)}.offsets(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        UnDirBoutDetails.BoutFirstSyllLabel(Index) = UnDirBoutDetails.BoutLabels{Index}(1);
        UnDirBoutDetails.SongBoutOrNot(Index) = SURecordingDetails.Bouts(j,7);
        
        Index = Index + 1;
    end
end

% Now do the same for directed song - only for directed song I will also
% consider all bouts that have enough time in the pre-bout window. I can
% always compare everything to spontaneous activity in undir song.

DirBouts = find((SURecordingDetails.BoutDirUnDir == 1));

if (isempty(DirBouts))
   DirBoutDetails = [];
end

% Now for each of the neurons, take the bouts with enough data at the
% beginning and find all the spikes before bout initiation in two different
% windows - -2 to -1.5s (spont. activity window) and -0.6 to -0.1s
% (pre-song window)
% In addition get the spike train, the raster and the pst for that bout.
% Keep track of bout # also

if (SURecordingDetails.Continuousdata == 1)
    Index = 1;

    for j = DirBouts(:)',
        
        DirBoutDetails.BoutIndex(Index) = j;
        
        DirBoutDetails.BoutLabels{Index} = AllLabels(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        DirBoutDetails.BoutOnsets{Index} = AllOnsets(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        DirBoutDetails.BoutOffsets{Index} = AllOffsets(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        DirBoutDetails.BoutFirstSyllLabel(Index) = DirBoutDetails.BoutLabels{Index}(1);
        DirBoutDetails.SongBoutOrNot(Index) = SURecordingDetails.Bouts(j,7);
        
        Index = Index + 1;
    end
else
    Index = 1;
    for j = DirBouts(:)',
        DirBoutDetails.BoutIndex(Index) = j;
        DirBoutDetails.BoutLabels{Index} = SURecordingDetails.NoteInfo{SURecordingDetails.Bouts(j,3)}.labels(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        DirBoutDetails.BoutOnsets{Index} = SURecordingDetails.NoteInfo{SURecordingDetails.Bouts(j,3)}.onsets(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        DirBoutDetails.BoutOffsets{Index} = SURecordingDetails.NoteInfo{SURecordingDetails.Bouts(j,3)}.offsets(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        
        DirBoutDetails.BoutFirstSyllLabel(Index) = DirBoutDetails.BoutLabels{Index}(1);
        
        DirBoutDetails.SongBoutOrNot(Index) = SURecordingDetails.Bouts(j,7);
        
        Index = Index + 1;
    end
end