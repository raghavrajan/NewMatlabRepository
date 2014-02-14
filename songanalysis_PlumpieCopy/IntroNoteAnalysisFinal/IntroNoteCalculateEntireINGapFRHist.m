function [INFR, GapFR, Edges, Position, Sig] = IntroNoteCalculateEntireINGapFRHist(Neural_INR, PreTime, PostTime, BinSize)

GapPreMotorLag = 0.045;

Width = 0.005;
GaussianLen = 2;
IFRFs = 1/BinSize;
XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (IFRFs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * IFRFs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * IFRFs) * (Width * IFRFs)));

Edges = -PreTime:BinSize:PostTime;

INIndex = 0;
for i = 1:length(Neural_INR.NoofINs),
    if (Neural_INR.NoofINs(i) > 0)
        BoutSpikeTimes = Neural_INR.BoutDetails(i).SpikeTimes;
        INs = Neural_INR.INs{i};
        for j = 1:length(INs),
            INIndex = INIndex + 1;
            Position(INIndex,1) = j - length(INs) - 1;
            Position(INIndex,2) = j;
            INOnset = Neural_INR.BoutDetails(i).onsets(INs(j));
            INOffset = Neural_INR.BoutDetails(i).offsets(INs(j));
            GapOnset = INOffset;
            GapOffset = Neural_INR.BoutDetails(i).onsets(INs(j)+1);
            
            INSpikeTimes{INIndex} = BoutSpikeTimes(find((BoutSpikeTimes >= (INOnset + Edges(1))) & (BoutSpikeTimes < (INOnset + Edges(end))))) - INOnset;
            GapSpikeTimes{INIndex} = BoutSpikeTimes(find((BoutSpikeTimes >= (GapOnset + Edges(1))) & (BoutSpikeTimes < (GapOnset + Edges(end))))) - INOnset;
            
            INFR(INIndex,:) = histc(BoutSpikeTimes, Edges + INOnset)/BinSize;
            GapFR(INIndex,:) = histc(BoutSpikeTimes, Edges + GapOnset)/BinSize;
        end
    end
end

for i = 1:size(Neural_INR.WithinBoutNoofINs, 1),
    if (Neural_INR.WithinBoutNoofINs(i,1) > 0)
        BoutSpikeTimes = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes;
        INs = Neural_INR.WithinBoutINs{i};
        for j = 1:length(INs),
            INIndex = INIndex + 1;
            Position(INIndex,1) = j - length(INs) - 1;
            Position(INIndex,2) = j;
            INOnset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(INs(j));
            INOffset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).offsets(INs(j));
            GapOnset = INOffset;
            GapOffset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(INs(j)+1);
                        
            INSpikeTimes{INIndex} = BoutSpikeTimes(find((BoutSpikeTimes >= (INOnset + Edges(1))) & (BoutSpikeTimes < (INOnset + Edges(end))))) - INOnset;
            GapSpikeTimes{INIndex} = BoutSpikeTimes(find((BoutSpikeTimes >= (GapOnset + Edges(1))) & (BoutSpikeTimes < (GapOnset + Edges(end))))) - INOnset;

            INFR(INIndex,:) = histc(BoutSpikeTimes, Edges + INOnset)/BinSize;
            GapFR(INIndex,:) = histc(BoutSpikeTimes, Edges + GapOnset)/BinSize;
        end
    end
end

INFR(:,end) = [];
GapFR(:,end) = [];
Edges = Edges(1:end-1) + BinSize/2;

HistogramDur = PostTime + PreTime;

ZeroTime = find(Edges > 0, 1, 'first');

for i = 1:1000,
    RandShift = rand(length(INSpikeTimes),1)*HistogramDur;
    PST = [];
    for j = 1:length(INSpikeTimes),
        SpikeTimes = INSpikeTimes{j} + RandShift(j);
        SpikeTimes(find(SpikeTimes >= PostTime)) = SpikeTimes(find(SpikeTimes >= PostTime)) - HistogramDur;
        PST(j,:) = histc(SpikeTimes, Edges)/BinSize;
    end
    MaxINPST(i) = max(mean(PST));
    MinINPST(i) = min(mean(PST));
end

for i = 1:1000,
    RandShift = rand(length(GapSpikeTimes),1)*HistogramDur;
    PST = [];
    for j = 1:length(GapSpikeTimes),
        SpikeTimes = GapSpikeTimes{j} + RandShift(j);
        SpikeTimes(find(SpikeTimes >= PostTime)) = SpikeTimes(find(SpikeTimes >= PostTime)) - HistogramDur;
        PST(j,:) = histc(SpikeTimes, Edges)/BinSize;
    end
    MaxGapPST(i) = max(mean(PST));
    MinGapPST(i) = min(mean(PST));
end

figure;
subplot(1,2,1);
plot(Edges, mean(INFR));
hold on;
plot([Edges(1) Edges(end)], [(mean(mean(INFR(:,1:round(0.1/BinSize))))) + 2*std(mean(INFR(:,1:round(0.1/BinSize)))) (mean(mean(INFR(:,1:round(0.1/BinSize))))) + 2*std(mean(INFR(:,1:round(0.1/BinSize))))], 'k--');plot([Edges(1) Edges(end)], [((mean(mean(INFR(:,1:round(0.1/BinSize))))) + 2*std(mean(INFR(:,1:round(0.1/BinSize))))) ((mean(mean(INFR(:,1:round(0.1/BinSize))))) + 2*std(mean(INFR(:,1:round(0.1/BinSize)))))], 'k--');
plot([Edges(1) Edges(end)], [(mean(mean(INFR(:,1:round(0.1/BinSize))))) - 2*std(mean(INFR(:,1:round(0.1/BinSize)))) (mean(mean(INFR(:,1:round(0.1/BinSize))))) - 2*std(mean(INFR(:,1:round(0.1/BinSize))))], 'k--');

subplot(1,2,2);
plot(Edges, mean(GapFR));
hold on;
plot([Edges(1) Edges(end)], [(mean(mean(GapFR(:,1:round(0.1/BinSize))))) + 2*std(mean(GapFR(:,1:round(0.1/BinSize)))) (mean(mean(GapFR(:,1:round(0.1/BinSize))))) + 2*std(mean(GapFR(:,1:round(0.1/BinSize))))], 'k--');
plot([Edges(1) Edges(end)], [(mean(mean(GapFR(:,1:round(0.1/BinSize))))) - 2*std(mean(GapFR(:,1:round(0.1/BinSize)))) (mean(mean(GapFR(:,1:round(0.1/BinSize))))) - 2*std(mean(GapFR(:,1:round(0.1/BinSize))))], 'k--');

MeanINPST = mean(INFR);
[MaxVal, MaxLoc] = max(MeanINPST(1:ZeroTime));
Sig(1,:) = [MaxVal Edges(MaxLoc) length(find(MaxINPST > MaxVal))/1000];

[MinVal, MinLoc] = min(MeanINPST(1:ZeroTime));
Sig(2,:) = [MinVal Edges(MinLoc) length(find(MinINPST < MinVal))/1000];

MeanGapPST = mean(GapFR);
[MaxVal, MaxLoc] = max(MeanGapPST(1:ZeroTime));
Sig(3,:) = [MaxVal Edges(MaxLoc) length(find(MaxGapPST > MaxVal))/1000];

[MinVal, MinLoc] = min(MeanGapPST(1:ZeroTime));
Sig(4,:) = [MinVal Edges(MinLoc) length(find(MinGapPST < MinVal))/1000];

disp('Finished');