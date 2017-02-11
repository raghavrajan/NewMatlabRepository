function [] = PlotSyllableSpikeOnsetOffset(FileList, NoteFileDir, SpikeFileDir, SpikeChan, SpikeCluster, Syllable)

Fid = fopen(FileList, 'r');
FileName = fgetl(Fid);

while (ischar(FileName(1)))
    ActualSpikeTimes = load([SpikeFileDir, '/', 'Chan', num2str(SpikeChan), '_', FileName, '.spk']);
    ActualSpikeTimes = ActualSpikeTimes(find(ActualSpikeTimes(:,1) == SpikeCluster), 2);
    NoteTimes = load([NoteFileDir, '/', FileName, '.not.mat']);
    
    Matches = find(NoteTimes.labels == Syllable);
        
    for i = 1:length(Matches),
        Range(1) = -0.25;
        Range(2) = 0.5;
        Resolution = 0.0001;
        
        Spikes = [];
        SpikeTrain = ActualSpikeTimes(find((ActualSpikeTimes >= (NoteTimes.onsets(Matches(i))/1000 - abs(Range(1)))) & (ActualSpikeTimes <= (NoteTimes.onsets(Matches(i))/1000 + Range(2)))));
        SpikeTrain = SpikeTrain - (NoteTimes.onsets(Matches(i))/1000);
        
              
        Time = Range(1):Resolution:Range(2);
        
        if ((SpikeTrain(1) ~= Range(1)) && (SpikeTrain(end) ~= Range(2)))
            SpikeTimes = [Range(1); [(SpikeTrain)]; Range(2)];
        else
            if (SpikeTrain(1) == Range(1))
                SpikeTimes = [[(SpikeTrain)]; Range(1)];
            else
                if (SpikeTrain(end) == Range(2))
                    SpikeTimes = [Range(1); [(SpikeTrain)]];
                else
                    SpikeTimes = [Range(1); [(SpikeTrain)]; Range(2)];
                end
            end
        end
        SpikeTimeIndices = round((SpikeTimes + abs(Range(1)))/Resolution);
        IFR = zeros(1,length(Time));
        for j = 2:length(SpikeTimeIndices),
            IFR(1,((SpikeTimeIndices(j-1) + 1):SpikeTimeIndices(j))) = 1/(abs(SpikeTimes(j) - SpikeTimes(j-1)));
        end
        plot(Time, IFR);
        hold on;
    end   
    
    FileName = fgetl(Fid);
end

fclose(Fid);