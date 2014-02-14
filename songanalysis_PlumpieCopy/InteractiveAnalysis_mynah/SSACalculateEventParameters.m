function [EventParameters] = SSACalculateEventParameters(FileInfo, MedianMotif, BinSize, PreSongStartDuration, PreSongEndDuration, FileInfo2, EventThreshold, EventWindowWidth, ContextString)

disp([ContextString, ':']);
EventWindowWidth = EventWindowWidth/2;

Edges = -PreSongStartDuration:BinSize:(MedianMotif.Length - PreSongEndDuration);

Flag = 1;

if (isfield(FileInfo2,'WEventParameters'))
    if (isfield(FileInfo2.WEventParameters, 'EventIndices'))
    EventParameters = FileInfo2.WEventParameters;
        if (isfield(FileInfo2.WEventParameters, 'EventTimes'))
            EventParameters = rmfield(EventParameters,'Jitter');
            EventParameters = rmfield(EventParameters, 'NoofSpikes');
            EventParameters = rmfield(EventParameters, 'ISIs');
            EventParameters = rmfield(EventParameters, 'FirstSpikeTime');
            EventParameters = rmfield(EventParameters, 'EventSpikeTimes');
            EventParameters = rmfield(EventParameters, 'STDFirstSpike');
            EventParameters = rmfield(EventParameters, 'JittertoMeanFirstSpikeTime');
            EventParameters = rmfield(EventParameters, 'RMSJitter');
            EventParameters = rmfield(EventParameters, 'EventTimes');
            EventParameters = rmfield(EventParameters, 'Reliability');
            EventParameters = rmfield(EventParameters, 'MeanValues');
        end
    else
        Flag = 0;
    end
else
    Flag = 0;
end

if (~Flag)
    EventParameters.EventThreshold = mean(mean(FileInfo.PST)) + (EventThreshold * std(mean(FileInfo.PST)));
    
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

    figure();
    for i = 1:length(FileInfo.WSpikeTrain),
        line([[FileInfo.WSpikeTrain{i}]'; [FileInfo.WSpikeTrain{i}]'], [(ones(size([FileInfo.WSpikeTrain{i}]))'*i + 0.25); (ones(size([FileInfo.WSpikeTrain{i}]))'*i - 0.25)], 'Color', 'k', 'LineWidth', 0.5);
    end
    axis tight;
    hold on;
    temp = axis;
    
    for i = 1:length(EventParameters.EventTimes),
        plot([(EventParameters.EventTimes(i) - EventWindowWidth) (EventParameters.EventTimes(i) - EventWindowWidth)], [temp(3) temp(4)], 'b');
        plot([(EventParameters.EventTimes(i) + EventWindowWidth) (EventParameters.EventTimes(i) + EventWindowWidth)], [temp(3) temp(4)], 'b');
    end
    axis([-PreSongStartDuration MedianMotif.Length temp(3) temp(4)]);
    
    for i = 1:length(FileInfo.WSpikeTrain),
        for j = 1:length(EventParameters.EventTimes),
            SpikeTrain = [FileInfo.WSpikeTrain{i}];
            UWSpikeTrain = [FileInfo.UWSpikeTrain{i}];
            EventSpikeIndices = find((SpikeTrain > (EventParameters.EventTimes(j) - EventWindowWidth)) & (SpikeTrain < (EventParameters.EventTimes(j) + EventWindowWidth)));
            if (length(EventSpikeIndices) > 0)
                EventSpikeTimes = SpikeTrain(EventSpikeIndices);
                UWEventSpikeTimes = UWSpikeTrain(EventSpikeIndices);
                if (length(UWEventSpikeTimes) > 1)
                    EventParameters.BurstDuration{i,j} = [diff(UWEventSpikeTimes) (ones(size(diff(UWEventSpikeTimes))) * (UWEventSpikeTimes(end) - UWEventSpikeTimes(1)))];
                else
                    EventParameters.BurstDuration{i,j} = 0;
                end
                EventParameters.Jitter(i,j) = EventSpikeTimes(1) - EventParameters.EventTimes(j);
                EventParameters.FirstSpikeTime(i,j) = EventSpikeTimes(1);
                EventParameters.UWFirstSpikeTime(i,j) = UWEventSpikeTimes(1);
                EventParameters.NoofSpikes(i,j) = length(EventSpikeTimes);
                EventParameters.ISIs{i,j} = diff(UWEventSpikeTimes);
                EventParameters.EventSpikeTimes{i,j} = EventSpikeTimes;
                EventParameters.Reliability(i,j) = 1;
            else
                EventParameters.Reliability(i,j) = 0;
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
    for i = 1:length(EventParameters.EventIndices),
        disp(['Event #',num2str(i),': ']);
        ISIs = [];
        for j = 1:size(EventParameters.ISIs,1),
            ISIs = [ISIs; EventParameters.ISIs{j,i}];
        end
        EventParameters.MeanValues.ISIs{i,1} = mean(ISIs);
        EventParameters.MeanValues.ISIs{i,2} = std(ISIs);
        EventParameters.MeanValues.Reliability{i} = mean(EventParameters.Reliability(:,i)) * 100;
        EventParameters.MeanValues.MeanNoofSpikes{i,1} = mean(EventParameters.NoofSpikes(:,i));
        EventParameters.MeanValues.MeanNoofSpikes{i,2} = std(EventParameters.NoofSpikes(:,i));
        disp(['Mean no of spikes is ', num2str(EventParameters.MeanValues.MeanNoofSpikes{i,1}), ' +/- ', num2str(EventParameters.MeanValues.MeanNoofSpikes{i,2})]);
        disp(['Reliabililty is ', num2str(EventParameters.MeanValues.Reliability{i}), '%']);
        disp(['Mean ISI is ', num2str(EventParameters.MeanValues.ISIs{i,1}), 's +/- ', num2str(EventParameters.MeanValues.ISIs{i,2}), 's']);
        disp(['RMS jitter to mean first spike time is ', num2str(EventParameters.JittertoMeanFirstSpikeTime(i)), 's']);
    end
end
