function [Correlation] = CalculateIFRCorr(SpikeTrain, Range, Resolution)

Time = Range(1):Resolution:Range(2);

IFR = zeros(length(SpikeTrain),length(Time));

for i = 1:length(SpikeTrain),
    SpikeTimeIndices = [0; round(([SpikeTrain{i}] + Range(1))/Resolution); length(Time)];
    if ((SpikeTrain{i}(1) ~= Range(1)) && (SpikeTrain{i}(end) ~= Range(2)))
        SpikeTimes = [Range(1); [(SpikeTrain{i})]; Range(2)];
    else
        if (SpikeTrain{i}(1) == Range(1))
            SpikeTimes = [[(SpikeTrain{i})]; Range(1)];
        else
            if (SpikeTrain{i}(end) == Range(2))
                SpikeTimes = [Range(1); [(SpikeTrain{i})]];
            else
                SpikeTimes = [Range(1); [(SpikeTrain{i})]; Range(2)];
            end
        end
    end
        
    for j = 2:length(SpikeTimeIndices),
        IFR(i,((SpikeTimeIndices(j-1) + 1):SpikeTimeIndices(j))) = 1/(abs(SpikeTimes(j) - SpikeTimes(j-1)));
    end
end

Correlation = 0;
for i = 1:size(IFR,1),
    for j = (i+1):size(IFR,1),
        if ((norm(IFR(i,:) - mean(IFR(i,:))) * norm(IFR(j,:) - mean(IFR(j,:)))) == 0)
            Temp(i,j) = 0;
        else
            Temp(i,j) = (((IFR(i,:) - mean(IFR(i,:))) * (IFR(j,:)' - mean(IFR(j,:))))/(norm(IFR(i,:) - mean(IFR(i,:))) * norm(IFR(j,:) - mean(IFR(j,:)))));
            Correlation = Correlation + (((IFR(i,:) - mean(IFR(i,:))) * (IFR(j,:)' - mean(IFR(j,:))))/(norm(IFR(i,:) - mean(IFR(i,:))) * norm(IFR(j,:) - mean(IFR(j,:)))));
        end
    end
end

Correlation = (Correlation * 2)/(size(IFR,1) * (size(IFR,1) - 1));

disp([ContextString, ': Resolution is ', num2str(Resolution), ' ms and the IFR Correlation = ',num2str(mean(Temp(find(Temp)))), ' +/- ', num2str(std(Temp(find(Temp))))]);
Correlation(1) = mean(Temp(find(Temp)));
Correlation(2) = std(Temp(find(Temp)));
