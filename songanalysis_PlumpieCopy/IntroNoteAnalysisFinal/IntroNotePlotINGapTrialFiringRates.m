function [] = IntroNotePlotINGapTrialFiringRates(Neural_INR)

PreTime = 0.05;
PostTime = 0;

IFRFs = 2000;

PreTimeIndex = round(PreTime*IFRFs);
PostTimeIndex = round(PostTime*IFRFs);

Width = 0.005;
GaussianLen = 2;
XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (IFRFs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * IFRFs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * IFRFs) * (Width * IFRFs)));

for i = 1:length(Neural_INR.BoutDetails),
    ST{i}(1,:) = conv(Neural_INR.BoutDetails(i).IFR(2,:), GaussWin);
end

MaxINs = max([max(Neural_INR.NoofINs) max(Neural_INR.WithinBoutNoofINs(:,1))]);

INIndex = zeros(1,MaxINs);

OnsetINSpikeTimesDur = [];
OnsetGapSpikeTimesDur = [];

OffsetINSpikeTimesDur = [];
OffsetGapSpikeTimesDur = [];

OnsetOffsetINDist = [];
OnsetOffsetGapDist = [];

for i = 1:length(Neural_INR.NoofINs),
    if (Neural_INR.NoofINs(i) > 0)
        INs = Neural_INR.INs{i};
        IFR = Neural_INR.BoutDetails(i).IFR;
        SpikeTimes = Neural_INR.BoutDetails(i).SpikeTimes;
        
        for j = 1:length(INs),
            INIndex(length(INs) - j + 1) = INIndex(length(INs) - j + 1) + 1;
            Index = INIndex(length(INs) - j + 1);
            
            INOnset = Neural_INR.BoutDetails(i).onsets(INs(j));
            INOffset = Neural_INR.BoutDetails(i).offsets(INs(j));
            GapOnset = INOffset;
            GapOffset = Neural_INR.BoutDetails(i).onsets(INs(j) + 1);
         
            INDur{length(INs) - j + 1}(Index) = INOffset - INOnset;
            GapDur{length(INs) - j + 1}(Index) = GapOffset - GapOnset;
            
            INSpikeTimes = SpikeTimes(find((SpikeTimes >= (INOnset - PreTime)) & (SpikeTimes < (INOffset + PostTime))));
            INSpikeTimes = INSpikeTimes(:);
            GapSpikeTimes = SpikeTimes(find((SpikeTimes >= (GapOnset - PreTime)) & (SpikeTimes < (GapOffset + PostTime))));
            GapSpikeTimes = GapSpikeTimes(:);
            
            OnsetOffsetINDist = [OnsetOffsetINDist; [(INSpikeTimes - INOnset)/(INOffset - INOnset) (INSpikeTimes - INOffset)/(INOffset - INOnset)]];
            OnsetOffsetGapDist = [OnsetOffsetGapDist; [(GapSpikeTimes - GapOnset)/(GapOffset - GapOnset) (GapSpikeTimes - GapOffset)/(GapOffset - GapOnset)]];
            
            OnsetINSpikeTimes = SpikeTimes(find((SpikeTimes >= (INOnset - PreTime)) & (SpikeTimes < (INOffset + PostTime)))) - INOnset;
            OnsetINSpikeTimes = OnsetINSpikeTimes(:);
            OnsetGapSpikeTimes = SpikeTimes(find((SpikeTimes >= (GapOnset - PreTime)) & (SpikeTimes < (GapOffset + PostTime)))) - GapOnset;
            OnsetGapSpikeTimes = OnsetGapSpikeTimes(:);
            
            OnsetINSpikeTimesDur = [OnsetINSpikeTimesDur; ones(size(OnsetINSpikeTimes))*(INOffset - INOnset) OnsetINSpikeTimes];
            OnsetGapSpikeTimesDur = [OnsetGapSpikeTimesDur; ones(size(OnsetGapSpikeTimes))*(GapOffset - GapOnset) OnsetGapSpikeTimes];
            
            OffsetINSpikeTimes = SpikeTimes(find((SpikeTimes >= (INOnset - PreTime)) & (SpikeTimes < (INOffset + PostTime)))) - INOffset;
            OffsetINSpikeTimes = OffsetINSpikeTimes(:);
            OffsetGapSpikeTimes = SpikeTimes(find((SpikeTimes >= (GapOnset - PreTime)) & (SpikeTimes < (GapOffset + PostTime)))) - GapOffset;
            OffsetGapSpikeTimes = OffsetGapSpikeTimes(:);
            
            OffsetINSpikeTimesDur = [OffsetINSpikeTimesDur; ones(size(OffsetINSpikeTimes))*(INOffset - INOnset) OffsetINSpikeTimes];
            OffsetGapSpikeTimesDur = [OffsetGapSpikeTimesDur; ones(size(OffsetGapSpikeTimes))*(GapOffset - GapOnset) OffsetGapSpikeTimes];
            
            INIFROnset = find(IFR(1,:) >= INOnset, 1, 'first');
            INIFROffset = find(IFR(1,:) >= INOffset, 1, 'first');
            
            GapIFROnset = INIFROffset;
            GapIFROffset = find(IFR(1,:) >= GapOffset, 1, 'first');
            
            TempINIFRTime = IFR(1,(INIFROnset - PreTimeIndex):(INIFROffset + PostTimeIndex)) - INOnset;
            TempINIFR = ST{i}(1, (INIFROnset - PreTimeIndex):(INIFROffset + PostTimeIndex));
            
            
            INIFR{length(INs) - j + 1}{Index} = [TempINIFRTime; TempINIFR];

            TempGapIFRTime = IFR(1,(GapIFROnset - PreTimeIndex):(GapIFROffset + PostTimeIndex)) - GapOnset;
            TempGapIFR = ST{i}(1, (GapIFROnset - PreTimeIndex):(GapIFROffset + PostTimeIndex));
            
            GapIFR{length(INs) - j + 1}{Index} = [TempGapIFRTime; TempGapIFR];
            
            TempGapIFRTime = IFR(1,(GapIFROnset - PreTimeIndex):(GapIFROffset + PostTimeIndex)) - GapOffset;
            TempGapIFR = ST{i}(1, (GapIFROnset - PreTimeIndex):(GapIFROffset + PostTimeIndex));
            
            OffsetGapIFR{length(INs) - j + 1}{Index} = [TempGapIFRTime; TempGapIFR];
            
            if (j == length(INs))
                FirstSyllOnset = GapOffset;
                FirstSyllOffset = Neural_INR.BoutDetails(i).offsets(INs(j) + 1);
                
                FirstSyllDur(Index) = FirstSyllOffset - FirstSyllOnset;
            
                FirstSyllIFROnset = find(IFR(1,:) >= FirstSyllOnset, 1, 'first');
                FirstSyllIFROffset = find(IFR(1,:) >= FirstSyllOffset, 1, 'first');
                
                TempFirstSyllIFRTime = IFR(1,(FirstSyllIFROnset - PreTimeIndex):(FirstSyllIFROffset + PostTimeIndex)) - FirstSyllOnset;
                TempFirstSyllIFR = ST{i}(1, (FirstSyllIFROnset - PreTimeIndex):(FirstSyllIFROffset + PostTimeIndex));
                
                FirstSyllIFR{Index} = [TempFirstSyllIFRTime; TempFirstSyllIFR];
            end
        end
    end
end

for i = 1:size(Neural_INR.WithinBoutNoofINs, 1),
    if (Neural_INR.WithinBoutNoofINs(i, 1) > 0)
        INs = Neural_INR.WithinBoutINs{i};
        IFR = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).IFR;
        SpikeTimes = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes;
        for j = 1:length(INs),
            INIndex(length(INs) - j + 1) = INIndex(length(INs) - j + 1) + 1;
            Index = INIndex(length(INs) - j + 1);
            
            INOnset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(INs(j));
            INOffset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).offsets(INs(j));
            GapOnset = INOffset;
            GapOffset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(INs(j) + 1);
            
            INDur{length(INs) - j + 1}(Index) = INOffset - INOnset;
            GapDur{length(INs) - j + 1}(Index) = GapOffset - GapOnset;
            
            INSpikeTimes = SpikeTimes(find((SpikeTimes >= (INOnset - PreTime)) & (SpikeTimes < (INOffset + PostTime))));
            INSpikeTimes = INSpikeTimes(:);
            GapSpikeTimes = SpikeTimes(find((SpikeTimes >= (GapOnset - PreTime)) & (SpikeTimes < (GapOffset + PostTime))));
            GapSpikeTimes = GapSpikeTimes(:);
            
            OnsetOffsetINDist = [OnsetOffsetINDist; [(INSpikeTimes - INOnset)/(INOffset - INOnset) (INSpikeTimes - INOffset)/(INOffset - INOnset)]];
            OnsetOffsetGapDist = [OnsetOffsetGapDist; [(GapSpikeTimes - GapOnset)/(GapOffset - GapOnset) (GapSpikeTimes - GapOffset)/(GapOffset - GapOnset)]];
            
            OnsetINSpikeTimes = SpikeTimes(find((SpikeTimes >= (INOnset - PreTime)) & (SpikeTimes < (INOffset + PostTime)))) - INOnset;
            OnsetINSpikeTimes = OnsetINSpikeTimes(:);
            OnsetGapSpikeTimes = SpikeTimes(find((SpikeTimes >= (GapOnset - PreTime)) & (SpikeTimes < (GapOffset + PostTime)))) - GapOnset;
            OnsetGapSpikeTimes = OnsetGapSpikeTimes(:);
            
            OnsetINSpikeTimesDur = [OnsetINSpikeTimesDur; ones(size(OnsetINSpikeTimes))*(INOffset - INOnset) OnsetINSpikeTimes];
            OnsetGapSpikeTimesDur = [OnsetGapSpikeTimesDur; ones(size(OnsetGapSpikeTimes))*(GapOffset - GapOnset) OnsetGapSpikeTimes];
            
            OffsetINSpikeTimes = SpikeTimes(find((SpikeTimes >= (INOnset - PreTime)) & (SpikeTimes < (INOffset + PostTime)))) - INOffset;
            OffsetINSpikeTimes = OffsetINSpikeTimes(:);
            OffsetGapSpikeTimes = SpikeTimes(find((SpikeTimes >= (GapOnset - PreTime)) & (SpikeTimes < (GapOffset + PostTime)))) - GapOffset;
            OffsetGapSpikeTimes = OffsetGapSpikeTimes(:);
            
            OffsetINSpikeTimesDur = [OffsetINSpikeTimesDur; ones(size(OffsetINSpikeTimes))*(INOffset - INOnset) OffsetINSpikeTimes];
            OffsetGapSpikeTimesDur = [OffsetGapSpikeTimesDur; ones(size(OffsetGapSpikeTimes))*(GapOffset - GapOnset) OffsetGapSpikeTimes];
            
            INIFROnset = find(IFR(1,:) >= INOnset, 1, 'first');
            INIFROffset = find(IFR(1,:) >= INOffset, 1, 'first');
            
            GapIFROnset = INIFROffset;
            GapIFROffset = find(IFR(1,:) >= GapOffset, 1, 'first');
            
            TempINIFRTime = IFR(1,(INIFROnset - PreTimeIndex):(INIFROffset + PostTimeIndex)) - INOnset;
            TempINIFR = ST{Neural_INR.WithinBoutINBoutIndices(i)}(1, (INIFROnset - PreTimeIndex):(INIFROffset + PostTimeIndex));
            
            INIFR{length(INs) - j + 1}{Index} = [TempINIFRTime; TempINIFR];
            
            TempGapIFRTime = IFR(1,(GapIFROnset - PreTimeIndex):(GapIFROffset + PostTimeIndex)) - GapOnset;
            TempGapIFR = ST{Neural_INR.WithinBoutINBoutIndices(i)}(1, (GapIFROnset - PreTimeIndex):(GapIFROffset + PostTimeIndex));
            
            GapIFR{length(INs) - j + 1}{Index} = [TempGapIFRTime; TempGapIFR];
            
            TempGapIFRTime = IFR(1,(GapIFROnset - PreTimeIndex):(GapIFROffset + PostTimeIndex)) - GapOffset;
            TempGapIFR = ST{Neural_INR.WithinBoutINBoutIndices(i)}(1, (GapIFROnset - PreTimeIndex):(GapIFROffset + PostTimeIndex));
            
            OffsetGapIFR{length(INs) - j + 1}{Index} = [TempGapIFRTime; TempGapIFR];
         
            if (j == length(INs))
                FirstSyllOnset = GapOffset;
                FirstSyllOffset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).offsets(INs(j) + 1);
                
                FirstSyllDur(Index) = FirstSyllOffset - FirstSyllOnset;
            
                FirstSyllIFROnset = find(IFR(1,:) >= FirstSyllOnset, 1, 'first');
                FirstSyllIFROffset = find(IFR(1,:) >= FirstSyllOffset, 1, 'first');
                
                TempFirstSyllIFRTime = IFR(1,(FirstSyllIFROnset - PreTimeIndex):(FirstSyllIFROffset + PostTimeIndex)) - FirstSyllOnset;
                TempFirstSyllIFR = ST{Neural_INR.WithinBoutINBoutIndices(i)}(1, (FirstSyllIFROnset - PreTimeIndex):(FirstSyllIFROffset + PostTimeIndex));
                
                FirstSyllIFR{Index} = [TempFirstSyllIFRTime; TempFirstSyllIFR];
            end
        end
    end
end

TrialNos = cellfun(@length, INIFR);

PlotNoofINs = min(5, length(find(TrialNos >= 5)));

figure;
Colours = ['rgbcm'];
AxisLims = [];
ColorAxis = [];
for i = 1:PlotNoofINs,
    subplot(PlotNoofINs, 3, ((i-1)*3 + 1));
    [Durs, Indices] = sort(cellfun(@length, INIFR{i}));
    Index = 0;
    MaxLen = size(INIFR{i}{Indices(end)}, 2);
    TempSortedINIFR = ones(length(Indices), MaxLen) * NaN;
    for j = Indices,
        Index = Index + 1;
        TempSortedINIFR(Index,1:size(INIFR{i}{j},2)) = INIFR{i}{j}(2,:);
        %plot(INIFR{i}{j}(1,:), INIFR{i}{j}(2,:)+500*(j-1), 'b');
        hold on;
    end
    contourf(INIFR{i}{j}(1,:), 1:1:length(Indices), TempSortedINIFR);
    axis tight;
    Temp = axis;
    AxisLims = [AxisLims; Temp];
    Temp = caxis;
    ColorAxis = [ColorAxis; Temp];
    
	subplot(PlotNoofINs, 3, ((i-1)*3 + 2):((i-1)*3 + 3));
    
    [Durs, Indices] = sort(cellfun(@length, GapIFR{i}));
    Index = 0;
    MaxLen = size(GapIFR{i}{Indices(end)}, 2);
    TempSortedGapIFR = ones(length(Indices), MaxLen) * NaN;
    for j = Indices,
        Index = Index + 1;
        TempSortedGapIFR(Index,1:size(GapIFR{i}{j},2)) = GapIFR{i}{j}(2,:);
        %plot(GapIFR{i}{j}(1,:), GapIFR{i}{j}(2,:)+500*(j-1), 'b');
        hold on;
    end
    contourf(GapIFR{i}{j}(1,:), 1:1:length(Indices), TempSortedGapIFR);
    axis tight;
    Temp = axis;
    AxisLims = [AxisLims; Temp];    
    Temp = caxis;
    ColorAxis = [ColorAxis; Temp];
end

for i = 1:PlotNoofINs,
    subplot(PlotNoofINs, 3, ((i-1)*3 + 1));
    axis([-PreTime max(AxisLims([1:2:size(AxisLims,1)], 2)) 0 AxisLims((i-1)*2 + 1, 4)]);
    caxis([0 max(ColorAxis(:,2))]);
    colorbar;
    
    subplot(PlotNoofINs, 3, ((i-1)*3 + 2):((i-1)*3 + 3));
    axis([-PreTime max(AxisLims([2:2:size(AxisLims,1)], 2)) 0 AxisLims((i-1)*2 + 2, 4)]);
    caxis([0 max(ColorAxis(:,2))]);
    colorbar;
end
% subplot(2, PlotNoofINs+1, PlotNoofINs + 1);
% [Durs, Indices] = sort(FirstSyllDur);
% for j = Indices,
%     plot(FirstSyllIFR{j}(1,:), FirstSyllIFR{j}(2,:), 'b');
%     hold on;
% end
% 
% axis tight;

figure;
Colours = ['rgbcm'];
AxisLims = [];
ColorAxis = [];
for i = 1:PlotNoofINs,
    subplot(PlotNoofINs, 2, ((i-1)*2 + 1));
    [Durs, Indices] = sort(cellfun(@length, GapIFR{i}));
    Index = 0;
    MaxLen = size(GapIFR{i}{Indices(end)}, 2);
    TempSortedGapIFR = ones(length(Indices), MaxLen) * NaN;
    for j = Indices,
        Index = Index + 1;
        TempSortedGapIFR(Index,1:size(GapIFR{i}{j},2)) = GapIFR{i}{j}(2,:);
        %plot(INIFR{i}{j}(1,:), INIFR{i}{j}(2,:)+500*(j-1), 'b');
        hold on;
    end
    contourf(GapIFR{i}{j}(1,:), 1:1:length(Indices), TempSortedGapIFR);
    axis tight;
    Temp = axis;
    AxisLims = [AxisLims; Temp];
    Temp = caxis;
    ColorAxis = [ColorAxis; Temp];
    
	subplot(PlotNoofINs, 2, ((i-1)*2 + 2));
    
    [Durs, Indices] = sort(cellfun(@length, OffsetGapIFR{i}));
    Index = 0;
    MaxLen = size(OffsetGapIFR{i}{Indices(end)}, 2);
    TempSortedOffsetGapIFR = ones(length(Indices), MaxLen) * NaN;
    for j = Indices,
        Index = Index + 1;
        TempSortedOffsetGapIFR(Index,MaxLen - size(OffsetGapIFR{i}{j},2) + 1:end) = OffsetGapIFR{i}{j}(2,:);
        %plot(GapIFR{i}{j}(1,:), GapIFR{i}{j}(2,:)+500*(j-1), 'b');
        hold on;
    end
    contourf(OffsetGapIFR{i}{j}(1,:), 1:1:length(Indices), TempSortedOffsetGapIFR);
    axis tight;
    Temp = axis;
    AxisLims = [AxisLims; Temp];    
    Temp = caxis;
    ColorAxis = [ColorAxis; Temp];
end

for i = 1:PlotNoofINs,
    subplot(PlotNoofINs, 2, ((i-1)*2 + 1));
    axis([-PreTime max(AxisLims([1:2:size(AxisLims,1)], 2)) 0 AxisLims((i-1)*2 + 1, 4)]);
    caxis([0 max(ColorAxis(:,2))]);
    colorbar;
    
    subplot(PlotNoofINs, 2, ((i-1)*2 + 2));
    axis([min(AxisLims([2:2:size(AxisLims,1)], 1)) 0 0 AxisLims((i-1)*2 + 2, 4)]);
    caxis([0 max(ColorAxis(:,2))]);
    colorbar;
end

disp('Finished');