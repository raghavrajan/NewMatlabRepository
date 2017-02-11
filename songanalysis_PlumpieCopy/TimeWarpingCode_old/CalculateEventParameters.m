function [EventParameters] = CalculateEventParameters(FileInfo, MedianMotif, BinSize, Latency, FileInfo2)

PreDataTime = 0.05;

Edges = -PreDataTime:BinSize:(MedianMotif.Length);

if (isfield(FileInfo2,'EventParameters'))
    EventParameters = FileInfo2.EventParameters;
    clear EventParameters.Jitter;
    clear EventParameters.NoofSpikes;
    clear EventParameters.ISIs;
else
    EventParameters.EventThreshold = mean(mean(FileInfo.PST)) + (4 * std(mean(FileInfo.PST)));
    
    TempEventIndices = find(mean(FileInfo.PST) > EventParameters.EventThreshold);
    if (length(TempEventIndices) > 0)
        EventIndices(:,1) = TempEventIndices;
        if ((size(EventIndices,1) == 1) && (size(EventIndices,2) > 1))
            EventIndices = EventIndices';
        end
    
        EventParameters.AllEventIndices = EventIndices;
        EventParameters.EventIndices = [EventIndices(1); EventIndices(find(diff(EventIndices) > 1) + 1)];
        EventParameters.ActualIndices = [0; (find(diff(EventIndices) > 1)); length(EventIndices)];
    end
end

if (isfield(EventParameters,'EventIndices'))
    for i = 1:length(EventParameters.EventIndices),
        EventParameters.EventTimes(i) = mean(Edges(EventParameters.AllEventIndices((EventParameters.ActualIndices(i) + 1):(EventParameters.ActualIndices(i+1)))));
    end

    for i = 1:length(FileInfo.SpikeTrain),
        for j = 1:length(EventParameters.EventTimes),
            SpikeTrain = [FileInfo.SpikeTrain{i} - Latency];
            EventSpikeIndices = find((SpikeTrain > (EventParameters.EventTimes(j) - 0.035)) & (SpikeTrain < (EventParameters.EventTimes(j) + 0.035)));
            if (length(EventSpikeIndices) > 0)
                EventSpikeTimes = SpikeTrain(EventSpikeIndices);
                EventParameters.Jitter(i,j) = EventSpikeTimes(1) - EventParameters.EventTimes(j);
                EventParameters.NoofSpikes(i,j) = length(EventSpikeTimes);
                EventParameters.ISIs{i,j} = diff(EventSpikeTimes);
                EventParameters.EventSpikeTimes{i,j} = EventSpikeTimes;
            end
        end
    end
    for i = 1:length(EventParameters.EventTimes),
        EventSpikeTimes = [];
        for j = 1:size(EventParameters.EventSpikeTimes,1),
            if (length(EventParameters.EventSpikeTimes{j,i}))
                EventSpikeTimes = [EventSpikeTimes [EventParameters.EventSpikeTimes{j,i}(1)]'];
            end
        end
        EventParameters.STDFirstSpike(i) = std(EventSpikeTimes);
        EventParameters.JittertoMeanFirstSpikeTime(i) = sqrt(sum((EventSpikeTimes - mean(EventSpikeTimes)).*(EventSpikeTimes - mean(EventSpikeTimes))))/sqrt(length(EventSpikeTimes));
        EventParameters.RMSJitter(i) = sqrt((sum(EventParameters.Jitter(find(EventParameters.Jitter(:,i)),i).*EventParameters.Jitter(find(EventParameters.Jitter(:,i)),i)))/(length(find(EventParameters.Jitter(:,i)))));                        
    end
end
