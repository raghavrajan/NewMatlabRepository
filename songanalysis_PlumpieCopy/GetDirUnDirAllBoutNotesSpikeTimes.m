function [BoutDetails] = GetDirUnDirAllBoutNotesSpikeTimes(SURecordingDetails)

% Get all bouts

SongBouts = 1:1:size(SURecordingDetails.Bouts(:,7),1); % I'm calling this song bout but it is essentially all bouts

% Now for each valid bout, get the onsets and offsets of motifs
% and the spike times within that bout

% Put all cluster spikes together into one long spike times list
SURecordingDetails.AllClusterSpikeTimes = [];
CumulativeFileTime = 0;
for i = 1:length(SURecordingDetails.ClusterSpikeTimes),
    SURecordingDetails.AllClusterSpikeTimes = [SURecordingDetails.AllClusterSpikeTimes; (SURecordingDetails.ClusterSpikeTimes{i}(:) + CumulativeFileTime)];
    CumulativeFileTime = CumulativeFileTime + SURecordingDetails.FileLen(i);
end

% Icnitialize dir undir parameters to NaN
BoutDetails.DirUnDir = ones(size(SURecordingDetails.Bouts,1),1)*NaN;

PrePostBoutTime = 100; % in ms
for j = SongBouts(:)',
    BoutDetails.DirUnDir(j) = SURecordingDetails.BoutDirUnDir(j);
    BoutDetails.SongBoutOrNot(j) = SURecordingDetails.Bouts(j,7);
    if (SURecordingDetails.Continuousdata == 1)
        BoutDetails.BoutLabels{j} = char(SURecordingDetails.SyllableData(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2), 1));
        BoutDetails.BoutOnsets{j} = SURecordingDetails.SyllableData(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2), 6);
        BoutDetails.BoutOffsets{j} = SURecordingDetails.SyllableData(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2), 7);
        BoutDetails.SpikeTimes{j} = SURecordingDetails.AllClusterSpikeTimes(find((SURecordingDetails.AllClusterSpikeTimes >= (BoutDetails.BoutOnsets{j}(1) - PrePostBoutTime)) & (SURecordingDetails.AllClusterSpikeTimes <= (BoutDetails.BoutOffsets{j}(end) + PrePostBoutTime))));
        BoutDetails.FileEndTime(j) = sum(SURecordingDetails.FileLen);
    else
        BoutDetails.BoutLabels{j} = SURecordingDetails.NoteInfo{SURecordingDetails.Bouts(j,3)}.labels(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        BoutDetails.BoutOnsets{j} = SURecordingDetails.NoteInfo{SURecordingDetails.Bouts(j,3)}.onsets(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        BoutDetails.BoutOffsets{j} = SURecordingDetails.NoteInfo{SURecordingDetails.Bouts(j,3)}.offsets(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        BoutDetails.SpikeTimes{j} = SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)}(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutDetails.BoutOnsets{j}(1) - PrePostBoutTime)) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} <= (BoutDetails.BoutOffsets{j}(end) + PrePostBoutTime))));
        BoutDetails.FileEndTime(j) = SURecordingDetails.FileLen(SURecordingDetails.Bouts(j,3));
    end
    BoutDetails.BoutSyllNos{j} = SURecordingDetails.BoutSyllNo(j,:);
    % Check whether there is enough pre-post time
    if ((BoutDetails.BoutOnsets{j}(1) < PrePostBoutTime) || (BoutDetails.BoutOffsets{j}(end) > (BoutDetails.FileEndTime(j) - PrePostBoutTime)))
        BoutDetails.EnoughData(j) = 0;
    else
        BoutDetails.EnoughData(j) = 1;
    end
end
