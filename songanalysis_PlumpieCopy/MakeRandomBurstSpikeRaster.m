function [] = MakeRandomBurstSpikeRaster(BurstProb)

Raster = [];
PST = [];
Edges = 0:1:1000;
TempSpikes = [];
for i = 1:100,
    for j = 1:50,
        Temp = rand;
        if (Temp < BurstProb(j))
            Temp = rand;
            %Temp = 0.15;
            Temp = Temp * 9;
            Temp = Temp + (j-1)*20;
            Raster((end + 1), :) = [Temp i];
            Raster((end + 1), :) = [Temp+5 i];
            Raster((end + 1), :) = [Temp+10 i];
            TempSpikes = [TempSpikes; Temp; Temp+5; Temp+10];
        end
    end
    PST(i,:) = histc(TempSpikes, Edges);
    PST(i,:) = conv(PST(i,:), hann(10), 'same');
    SpikeTrain{i} = TempSpikes/1000;
    TempSpikes = [];
end
PST = PST/0.001;

figure;
%errorbar(Edges, mean(PST),std(PST), 'k');
plot(Edges, mean(PST), 'k');
hold on;
Raster(:,2) = Raster(:,2) + max(std(PST)) +max(mean(PST));
PlotRaster(Raster);
axis tight;
MedianMotif.Length = 1;
[Corr] = SSACalculateCorrGaussSameSize(SpikeTrain, MedianMotif, 0.01, 0, 0, 'Directed', 4);
