function [Raster] = MakeSpikeBurstRaster(SpikeProb, BurstProb, NoofTrials, BinSize, SongLength, CorrYesNo)

Raster = [];
Edges = 0:BinSize:SongLength;
PST = [];
SpikeTrain = [];
XGauss = 1:1:(1 + round(2 * 1 * 5 * BinSize));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((5 * BinSize) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (5 * BinSize) * (5 * BinSize)));

for i = 1:NoofTrials,
    TempSpikes = [];
    j = 0;
    while (j <= SongLength),
        SpikeOrNoSpike = rand;
        if (SpikeOrNoSpike < SpikeProb)
            SpikeOrBurst = rand;
            for k = 1:length(BurstProb),
                if (SpikeOrBurst <= BurstProb(k))
                    NoofSpikes = k;
                    break;
                end
            end
            for k = 1:NoofSpikes,
                ISI = rand;
                ISI = ISI*2 - 1;
                Raster(end+1,:) = [(j + ((k-1)*4 + ISI)) i];
                if (Raster(end, 1) < 0)
                    Raster(end, 1) = 0;
                end
                
                TempSpikes = [TempSpikes; Raster(end,1)];
            end
            if (k > 1)
                j = Raster(end,1) + 5;                
            else
                j = Raster(end,1) + 2;
            end
        end
        j = j + 0.1;
    end
    if (length(TempSpikes) == 0)
        PST(i,:) = zeros(1, length(Edges));
    else
        PST(i,:) = histc(TempSpikes, Edges);
    end
    CPST(i,:) = conv(PST(i,:), [1 1 1 1 1 1 1 1 1 1], 'same');
    HPST(i,:) = conv(PST(i,:), hann(10), 'same');
    GPST(i,:) = conv(PST(i,:), GaussWin, 'same');
    SpikeTrain{i} = TempSpikes/1000;
end
%CPST = PST;
CPST = CPST/(BinSize/1000);
CPST = CPST;
HPST = HPST/(BinSize/1000);
HPST = HPST;
GPST = GPST/(BinSize/1000);
GPST = GPST;
Edges(end) = [];
CPST(:,end) = [];
PST(:,end) = [];
HPST(:,end) = [];
GPST(:,end) = [];
figure;
subplot(2,1,2);
%plot(Edges, mean(CPST), 'k');
hold on;
%plot(Edges, mean(HPST), 'r');
plot(Edges, mean(GPST), 'b');

%axis([0 SongLength 0 max(max(mean(CPST)), max(mean(HPST)))]);
axis([0 SongLength 0 max(max(mean(GPST)), max(mean(GPST)))]);
%Raster(:,2) = Raster(:,2) + max(mean(CPST));
subplot(2,1,1);
line([Raster(:,1) Raster(:,1)]', [(Raster(:,2) + 0.4) (Raster(:,2) - 0.4)]', 'Color', 'k', 'MarkerSize', 6)
axis([0 SongLength 0 (Raster(end,2) + 1)]);
p = anova1(PST, [], 'off');
disp(['P value is ', num2str(p)]);
MedianMotif.Length = SongLength/1000;
if (strfind(CorrYesNo, 'yes'))
    Corr = [];
    for i = [1 5 10 25 50],
        [Corr(end+1,:)] = SSACalculateCorrGaussSameSize(SpikeTrain, MedianMotif, i/1000, 0, 0, 'Simulated', 4);
    end
    figure
    plot(Corr(:,1)*1000, Corr(:,2), 'ks-');
end

FileInfo.UWSpikeTrain = SpikeTrain;
FileInfo.SongLengths = ones(length(SpikeTrain), 1) * SongLength/1000;
[MotifFiringRate, BurstFraction] = SSACalculateMotifFiringRate(FileInfo, 0, 0, 'Simulated');
