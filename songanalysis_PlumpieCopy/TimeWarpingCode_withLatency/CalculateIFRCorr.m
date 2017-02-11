function [Correlation] = CalculateIFRCorr(SpikeTrain, MedianMotif, Latency)

%Latency = 0.1;

Time = 0:0.0001:(MedianMotif.Length + Latency);

IFR = zeros(length(SpikeTrain),length(Time));

for i = 1:length(SpikeTrain),
    SpikeTimeIndices = [0; round(([SpikeTrain{i}] + Latency)/0.0001); length(Time)];
    SpikeTimes = [0; [(SpikeTrain{i} + Latency)]; (MedianMotif.Length + Latency)];
    for j = 2:length(SpikeTimeIndices),
        IFR(i,((SpikeTimeIndices(j-1) + 1):SpikeTimeIndices(j))) = 1/(SpikeTimes(j) - SpikeTimes(j-1));
    end
end

Correlation = 0;
for i = 1:size(IFR,1),
    for j = (i+1):size(IFR,1),
        Correlation = Correlation + (((IFR(i,:) - mean(IFR(i,:))) * (IFR(j,:)' - mean(IFR(j,:))))/(norm(IFR(i,:) - mean(IFR(i,:))) * norm(IFR(j,:) - mean(IFR(j,:)))));
    end
end

Correlation = (Correlation * 2)/(size(IFR,1) * (size(IFR,1) - 1));
disp(['Correlation = ',num2str(Correlation)]);