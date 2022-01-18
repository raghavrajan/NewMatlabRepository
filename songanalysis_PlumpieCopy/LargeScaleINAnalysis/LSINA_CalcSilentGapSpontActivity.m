function [UnDirGapSpontActivity, DirGapSpontActivity] = LSINA_CalcSilentGapSpontActivity(SURecordingDetails)

SpontActivityWindowSize = 500; % in milliseconds
SilentPeriodPrePostWindowSize = 3000; % in milliseconds

% Now if the data is continuous, then put it together into one long stream
% of spikes and one long stream of syllables

if (SURecordingDetails.Continuousdata == 1)
    ClusterSpikeData = [];
    AllClusterSpikeData = []; % includes optional cluster times as well
    TotalFileLen = 0;

    % First get all spike times into one long array
    for j = 1:length(SURecordingDetails.ClusterSpikeTimes),
        if (~isempty(SURecordingDetails.ClusterSpikeTimes{j}))
            ClusterSpikeData = [ClusterSpikeData; (TotalFileLen + SURecordingDetails.ClusterSpikeTimes{j}(:))];
        end
        if (~isempty(SURecordingDetails.OptionalClusterSpikeTimes{j}))
            AllClusterSpikeData = [AllClusterSpikeData; (TotalFileLen + SURecordingDetails.OptionalClusterSpikeTimes{j}(:))];
        end
        TotalFileLen = TotalFileLen + SURecordingDetails.FileLen(j);
    end

    % Now I want to take all 0.5s periods that have atleast 3s of
    % silence before and after. So, basically look for any silent
    % period that is longer than 6.5s and then pull out all the
    % non-overlapping 0.5s periods of spontaneous activity that I can
    % find within these. 
    % I will pull these out from the middle of the silent period. So
    % for eg: if there is 7.7s of silence, then I will pull out 3 0.5s
    % periods from the middle of this as follows - 3.1s silence - Gap1
    % (0.5s) - Gap2 (0.5s) - Gap3 (0.5s) - 3.1s silence.

    % Now, will only consider # of spikes in each 0.5s period and store
    % them in one long array.

    SilentGaps = find(SURecordingDetails.Gaps(:,1) >= (2*SilentPeriodPrePostWindowSize + SpontActivityWindowSize));
    SURecordingDetails.UnDirClusterSpontPeriodSpikeCount = [];
    SURecordingDetails.UnDirAllClusterSpontPeriodSpikeCount = [];
    SURecordingDetails.DirClusterSpontPeriodSpikeCount = [];
    SURecordingDetails.DirAllClusterSpontPeriodSpikeCount = [];
    for j = 1:length(SilentGaps),
        NumSilentPeriods = floor((SURecordingDetails.Gaps(SilentGaps(j),1) - 2*SilentPeriodPrePostWindowSize)/SpontActivityWindowSize);
        ExtraSilentPeriod = SURecordingDetails.Gaps(SilentGaps(j),1) - 2*SilentPeriodPrePostWindowSize - NumSilentPeriods*SpontActivityWindowSize;
        SilentPeriodStart = SURecordingDetails.Gaps(SilentGaps(j),6) + SilentPeriodPrePostWindowSize + ExtraSilentPeriod/2;
        for k = 1:NumSilentPeriods,
            % CHeck to see if this silent gap is within directed
            % presentation or not
            switch (SURecordingDetails.GapDirUnDir(SilentGaps(j)))
                case 1
                    SURecordingDetails.DirClusterSpontPeriodSpikeCount(end+1) = length(find((ClusterSpikeData >= (SilentPeriodStart + ((k-1)*SpontActivityWindowSize))) & (ClusterSpikeData < (SilentPeriodStart + (k*SpontActivityWindowSize)))));
                    SURecordingDetails.DirAllClusterSpontPeriodSpikeCount(end+1) = SURecordingDetails.DirClusterSpontPeriodSpikeCount(end) + length(find((AllClusterSpikeData >= (SilentPeriodStart + ((k-1)*SpontActivityWindowSize))) & (AllClusterSpikeData < (SilentPeriodStart + (k*SpontActivityWindowSize)))));

                case 0
                    SURecordingDetails.UnDirClusterSpontPeriodSpikeCount(end+1) = length(find((ClusterSpikeData >= (SilentPeriodStart + ((k-1)*SpontActivityWindowSize))) & (ClusterSpikeData < (SilentPeriodStart + (k*SpontActivityWindowSize)))));
                    SURecordingDetails.UnDirAllClusterSpontPeriodSpikeCount(end+1) = SURecordingDetails.UnDirClusterSpontPeriodSpikeCount(end) + length(find((AllClusterSpikeData >= (SilentPeriodStart + ((k-1)*SpontActivityWindowSize))) & (AllClusterSpikeData < (SilentPeriodStart + (k*SpontActivityWindowSize)))));
            end

        end
    end
else
    % for triggered data, I can look at all gaps and see if I manage to
    % find gaps that are longer than the gap length specified above and
    % get the spontaneous activity in such gaps
    SilentGaps = find(SURecordingDetails.Gaps(:,1) >= (2*SilentPeriodPrePostWindowSize + SpontActivityWindowSize));
    SURecordingDetails.UnDirClusterSpontPeriodSpikeCount = [];
    SURecordingDetails.UnDirAllClusterSpontPeriodSpikeCount = [];
    SURecordingDetails.DirClusterSpontPeriodSpikeCount = [];
    SURecordingDetails.DirAllClusterSpontPeriodSpikeCount = [];

    for j = 1:length(SilentGaps),
        NumSilentPeriods = floor((SURecordingDetails.Gaps(SilentGaps(j),1) - 2*SilentPeriodPrePostWindowSize)/SpontActivityWindowSize);
        ExtraSilentPeriod = SURecordingDetails.Gaps(SilentGaps(j),1) - 2*SilentPeriodPrePostWindowSize - NumSilentPeriods*SpontActivityWindowSize;
        SilentPeriodStart = SURecordingDetails.Gaps(SilentGaps(j),2) + SilentPeriodPrePostWindowSize + ExtraSilentPeriod/2;
        for k = 1:NumSilentPeriods,
            % CHeck to see if this silent gap is within directed
            % presentation or not
            switch (SURecordingDetails.GapDirUnDir(SilentGaps(j)))
                case 1
                    SURecordingDetails.DirClusterSpontPeriodSpikeCount(end+1) = length(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Gaps(SilentGaps(j),4)} >= (SilentPeriodStart + ((k-1)*SpontActivityWindowSize))) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Gaps(SilentGaps(j),4)} < (SilentPeriodStart + (k*SpontActivityWindowSize)))));
                    SURecordingDetails.DirAllClusterSpontPeriodSpikeCount(end+1) = SURecordingDetails.DirClusterSpontPeriodSpikeCount(end) + length(find((SURecordingDetails.OptionalClusterSpikeTimes{SURecordingDetails.Gaps(SilentGaps(j),4)} >= (SilentPeriodStart + ((k-1)*SpontActivityWindowSize))) & (SURecordingDetails.OptionalClusterSpikeTimes{SURecordingDetails.Gaps(SilentGaps(j),4)} < (SilentPeriodStart + (k*SpontActivityWindowSize)))));

                case 0
                    SURecordingDetails.UnDirClusterSpontPeriodSpikeCount(end+1) = length(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Gaps(SilentGaps(j),4)} >= (SilentPeriodStart + ((k-1)*SpontActivityWindowSize))) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Gaps(SilentGaps(j),4)} < (SilentPeriodStart + (k*SpontActivityWindowSize)))));
                    SURecordingDetails.UnDirAllClusterSpontPeriodSpikeCount(end+1) = SURecordingDetails.UnDirClusterSpontPeriodSpikeCount(end) + length(find((SURecordingDetails.OptionalClusterSpikeTimes{SURecordingDetails.Gaps(SilentGaps(j),4)} >= (SilentPeriodStart + ((k-1)*SpontActivityWindowSize))) & (SURecordingDetails.OptionalClusterSpikeTimes{SURecordingDetails.Gaps(SilentGaps(j),4)} < (SilentPeriodStart + (k*SpontActivityWindowSize)))));
            end
        end
    end
end

UnDirGapSpontActivity = SURecordingDetails.UnDirClusterSpontPeriodSpikeCount;
DirGapSpontActivity = SURecordingDetails.DirClusterSpontPeriodSpikeCount;
