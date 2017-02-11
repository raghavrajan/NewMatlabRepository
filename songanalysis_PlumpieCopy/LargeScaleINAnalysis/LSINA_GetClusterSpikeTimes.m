function [ClusterSpikeTimes, OptionalClusterSpikeTimes] = LSINA_GetClusterSpikeTimes(SURecordingDetails)

for i = 1:length(SURecordingDetails.SpikeInfo),
    ClusterSpikeTimes{i} = [];
    OptionalClusterSpikeTimes{i} = [];
    for j = 1:length(SURecordingDetails.SpikeClusterNos),
        ClusterSpikeIndices = find(SURecordingDetails.SpikeInfo{i}(:,1) == SURecordingDetails.SpikeClusterNos(j));
        ClusterSpikeTimes{i} = [ClusterSpikeTimes{i}; SURecordingDetails.SpikeInfo{i}(ClusterSpikeIndices(:),2)];
    end
    
    for j = 1:length(SURecordingDetails.Optionalclusterstotest),
        OptionalClusterSpikeIndices = find(SURecordingDetails.SpikeInfo{i}(:,1) == SURecordingDetails.Optionalclusterstotest(j));
        OptionalClusterSpikeTimes{i} = [OptionalClusterSpikeTimes{i}; SURecordingDetails.SpikeInfo{i}(OptionalClusterSpikeIndices(:),2)];
    end
end