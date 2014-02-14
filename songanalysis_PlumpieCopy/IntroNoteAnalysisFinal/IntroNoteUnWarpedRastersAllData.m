function [TempAllFR, TempAllFirstINCorr, TempAllLastINCorr, TempAllFFWinValue, LastvsOtherFRDiff, OthervsOtherFRDiff, LastvsOtherCorrDiff, OthervsOtherCorrDiff, Sig] = IntroNoteUnWarpedRastersAllData(Neural_INR, BinSize, MinNumber, MinPlotNumber, PreTime, PostTime, varargin)

LastvsOtherCorrDiff = [];
LastvsOtherFRDiff = [];
OthervsOtherCorrDiff = [];
OthervsOtherFRDiff = [];
            
if (nargin > 6)
    if (~isempty(varargin{1}))
        PlotOption = varargin{1};
    else
        PlotOption = 'on';
    end
    if (nargin > 7)
        ComparisonIN = varargin{2};
    else
        ComparisonIN = -1;
    end
else
    PlotOption = 'on';
    ComparisonIN = -1;
end

FFWindow = 0.03;
FFAdvance = 0.001;

Width = 0.001;
GaussianLen = 2;
IFRFs = 1/BinSize;
XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (IFRFs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * IFRFs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * IFRFs) * (Width * IFRFs)));

Edges = -PreTime:BinSize:PostTime;
FFEdges = Edges(1):FFAdvance:Edges(end);

FFDataStructEdges = (FFEdges(1) - FFWindow):0.001:(FFEdges(end)+FFWindow);

GapPreMotorLag = 0.045;

SyllDur = [];
INIndex = 0;
for i = 1:length(Neural_INR.NoofINs),
    if (Neural_INR.NoofINs(i) > 0)
        BoutSpikeTimes = Neural_INR.BoutDetails(i).SpikeTimes;
        INs = Neural_INR.INs{i};
        for j = 1:length(INs),
            INIndex = INIndex + 1;
            
            Position(INIndex,1) = j - length(INs) - 1;
            Position(INIndex,2) = j;
            NoofINs(INIndex,1) = length(INs);
            
            INOnset = Neural_INR.BoutDetails(i).onsets(INs(j));
            INOffset = Neural_INR.BoutDetails(i).offsets(INs(j));
            GapOnset = INOffset;
            GapOffset = Neural_INR.BoutDetails(i).onsets(INs(j)+1);
            
            RasterSpikeIndices = find((BoutSpikeTimes >= (INOnset - PreTime)) & (BoutSpikeTimes < (INOnset + PostTime)));
            TempSpikeTimes = BoutSpikeTimes(RasterSpikeIndices) - INOnset;
            TempSpikeTimes = TempSpikeTimes(:);
            RasterSpikeTimes{INIndex} = TempSpikeTimes;
            PST(INIndex,:) = conv(histc(TempSpikeTimes, Edges)/BinSize, GaussWin, 'same');
            
            INGapSpikeIndices = find((BoutSpikeTimes >= (INOnset - GapPreMotorLag)) & (BoutSpikeTimes < (GapOffset - GapPreMotorLag)));
            INGapFR(INIndex) = length(INGapSpikeIndices)/(GapOffset - INOnset);
            INGapNumSpike(INIndex) = length(INGapSpikeIndices);
            
            INDur(INIndex) = INOffset - INOnset;
            GapDur(INIndex) = GapOffset - GapOnset;
            INGapDur(INIndex) = GapOffset - INOnset;
            
            if (j == length(INs))
                SyllOnset = Neural_INR.BoutDetails(i).onsets(INs(j)+1);
                SyllOffset = Neural_INR.BoutDetails(i).offsets(INs(j)+1);
                SyllDur = [SyllDur; (SyllOffset - SyllOnset)];
            end
            
            FFSpikeIndices = find((BoutSpikeTimes >= (INOnset + FFEdges(1) - FFWindow)) & (BoutSpikeTimes < (INOnset + FFEdges(end) + FFWindow)));
            TempSpikeTimes = BoutSpikeTimes(FFSpikeIndices) - INOnset;
            TempSpikeTimes = TempSpikeTimes(:);
            FFDataStruct(INIndex,:) = histc(TempSpikeTimes, FFDataStructEdges);
            
            FFSpikeTimes{INIndex} = TempSpikeTimes;
            for FFWinIndex = 1:length(FFEdges),
               FFPST(INIndex, FFWinIndex) = length(find((FFSpikeTimes{INIndex} >= (FFEdges(FFWinIndex) - FFWindow/2)) & (FFSpikeTimes{INIndex} < (FFEdges(FFWinIndex) + FFWindow/2))));
            end
        end
        
        j = length(INs) + 1;
        INIndex = INIndex + 1;
            
        Position(INIndex,1) = j - length(INs) - 1;
        Position(INIndex,2) = j;
        NoofINs(INIndex,1) = length(INs);

        INOnset = Neural_INR.BoutDetails(i).onsets(INs(end) + 1);
        INOffset = Neural_INR.BoutDetails(i).offsets(INs(end) + 1);
        GapOnset = INOffset;
        if (length(Neural_INR.BoutDetails(i).onsets) >= (INs(end)+2))
            GapOffset = Neural_INR.BoutDetails(i).onsets(INs(end)+2);
        else
            GapOffset = GapOnset + 2;
        end
        
        RasterSpikeIndices = find((BoutSpikeTimes >= (INOnset - PreTime)) & (BoutSpikeTimes < (INOnset + PostTime)));
        TempSpikeTimes = BoutSpikeTimes(RasterSpikeIndices) - INOnset;
        TempSpikeTimes = TempSpikeTimes(:);
        RasterSpikeTimes{INIndex} = TempSpikeTimes;
        PST(INIndex,:) = conv(histc(TempSpikeTimes, Edges)/BinSize, GaussWin, 'same');

        INGapSpikeIndices = find((BoutSpikeTimes >= (INOnset - GapPreMotorLag)) & (BoutSpikeTimes < (GapOffset - GapPreMotorLag)));
        INGapFR(INIndex) = length(INGapSpikeIndices)/(GapOffset - INOnset);
        INGapNumSpike(INIndex) = length(INGapSpikeIndices);

        INSpikeIndices = find((BoutSpikeTimes >= (INOnset - GapPreMotorLag)) & (BoutSpikeTimes < (INOffset - GapPreMotorLag)));
        INFR(INIndex) = length(INSpikeIndices)/(INOffset - INOnset);
        INNumSpike(INIndex) = length(INSpikeIndices);
        
        GapSpikeIndices = find((BoutSpikeTimes >= (GapOnset - GapPreMotorLag)) & (BoutSpikeTimes < (GapOffset - GapPreMotorLag)));
        GapFR(INIndex) = length(GapSpikeIndices)/(GapOffset - GapOnset);
        GapNumSpike(INIndex) = length(GapSpikeIndices);
        
        INDur(INIndex) = INOffset - INOnset;
        INGapDur(INIndex) = GapOffset - INOnset;

        FFSpikeIndices = find((BoutSpikeTimes >= (INOnset + FFEdges(1) - FFWindow)) & (BoutSpikeTimes < (INOnset + FFEdges(end) + FFWindow)));
        TempSpikeTimes = BoutSpikeTimes(FFSpikeIndices) - INOnset;
        TempSpikeTimes = TempSpikeTimes(:);
        FFDataStruct(INIndex,:) = histc(TempSpikeTimes, FFDataStructEdges);

        FFSpikeTimes{INIndex} = TempSpikeTimes;
        for FFWinIndex = 1:length(FFEdges),
           FFPST(INIndex, FFWinIndex) = length(find((FFSpikeTimes{INIndex} >= (FFEdges(FFWinIndex) - FFWindow/2)) & (FFSpikeTimes{INIndex} < (FFEdges(FFWinIndex) + FFWindow/2))));
        end
      
    end
end

for i = 1:size(Neural_INR.WithinBoutNoofINs, 1),
    if (Neural_INR.WithinBoutNoofINs(i,1) > 0)
        BoutSpikeTimes = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes;
        INs = Neural_INR.WithinBoutINs{i};
        for j = 1:length(INs),
            INIndex = INIndex + 1;
            INSpikeIndices = find((BoutSpikeTimes >= (INOnset - GapPreMotorLag)) & (BoutSpikeTimes < (INOffset - GapPreMotorLag)));
            INFR(INIndex) = length(INSpikeIndices)/(INOffset - INOnset);
            INNumSpike(INIndex) = length(INSpikeIndices);

            GapSpikeIndices = find((BoutSpikeTimes >= (GapOnset - GapPreMotorLag)) & (BoutSpikeTimes < (GapOffset - GapPreMotorLag)));
            GapFR(INIndex) = length(GapSpikeIndices)/(GapOffset - GapOnset);
            GapNumSpike(INIndex) = length(GapSpikeIndices);Position(INIndex,1) = j - length(INs) - 1;
            Position(INIndex,2) = j;
            NoofINs(INIndex,1) = length(INs);
            
            INOnset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(INs(j));
            INOffset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).offsets(INs(j));
            GapOnset = INOffset;
            GapOffset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(INs(j)+1);
            
            RasterSpikeIndices = find((BoutSpikeTimes >= (INOnset - PreTime)) & (BoutSpikeTimes < (INOnset + PostTime)));
            TempSpikeTimes = BoutSpikeTimes(RasterSpikeIndices) - INOnset;
            TempSpikeTimes = TempSpikeTimes(:);
            RasterSpikeTimes{INIndex} = TempSpikeTimes;
            PST(INIndex,:) = conv(histc(TempSpikeTimes, Edges)/BinSize, GaussWin, 'same');
            
            INGapSpikeIndices = find((BoutSpikeTimes >= (INOnset - GapPreMotorLag)) & (BoutSpikeTimes < (GapOffset - GapPreMotorLag)));
            INGapFR(INIndex) = length(INGapSpikeIndices)/(GapOffset - INOnset);
            INGapNumSpike(INIndex) = length(INGapSpikeIndices);
            
            INDur(INIndex) = INOffset - INOnset;
            GapDur(INIndex) = GapOffset - INOnset;
            INGapDur(INIndex) = GapOffset - INOnset;
            
            if (j == length(INs))
                SyllOnset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(INs(j)+1);
                SyllOffset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).offsets(INs(j)+1);
                SyllDur = [SyllDur; (SyllOffset - SyllOnset)];
            end
            
            FFSpikeIndices = find((BoutSpikeTimes >= (INOnset + FFEdges(1) - FFWindow)) & (BoutSpikeTimes < (INOnset + FFEdges(end) + FFWindow)));
            TempSpikeTimes = BoutSpikeTimes(FFSpikeIndices) - INOnset;
            TempSpikeTimes = TempSpikeTimes(:);
            FFDataStruct(INIndex,:) = histc(TempSpikeTimes, FFDataStructEdges);
            
            FFSpikeTimes{INIndex} = TempSpikeTimes;
            for FFWinIndex = 1:length(FFEdges),
               FFPST(INIndex, FFWinIndex) = length(find((FFSpikeTimes{INIndex} >= (FFEdges(FFWinIndex) - FFWindow/2)) & (FFSpikeTimes{INIndex} < (FFEdges(FFWinIndex) + FFWindow/2))));
            end
        end
        
        j = length(INs) + 1;

        INIndex = INIndex + 1;
        Position(INIndex,1) = j - length(INs) - 1;
        Position(INIndex,2) = j;
        NoofINs(INIndex,1) = length(INs);

        INOnset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(INs(end) + 1);
        INOffset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).offsets(INs(end) + 1);
        GapOnset = INOffset;
        if (length(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets >= (INs(end) + 2)))
            GapOffset = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(INs(end)+2);
        else
            GapOffset = GapOnset + 0.2;
        end
        RasterSpikeIndices = find((BoutSpikeTimes >= (INOnset - PreTime)) & (BoutSpikeTimes < (INOnset + PostTime)));
        TempSpikeTimes = BoutSpikeTimes(RasterSpikeIndices) - INOnset;
        TempSpikeTimes = TempSpikeTimes(:);
        RasterSpikeTimes{INIndex} = TempSpikeTimes;
        PST(INIndex,:) = conv(histc(TempSpikeTimes, Edges)/BinSize, GaussWin, 'same');

        INGapSpikeIndices = find((BoutSpikeTimes >= (INOnset - GapPreMotorLag)) & (BoutSpikeTimes < (GapOffset - GapPreMotorLag)));
        INGapFR(INIndex) = length(INGapSpikeIndices)/(GapOffset - INOnset);
        INGapNumSpike(INIndex) = length(INGapSpikeIndices);

        INSpikeIndices = find((BoutSpikeTimes >= (INOnset - GapPreMotorLag)) & (BoutSpikeTimes < (INOffset - GapPreMotorLag)));
        INFR(INIndex) = length(INSpikeIndices)/(INOffset - INOnset);
        INNumSpike(INIndex) = length(INSpikeIndices);
        
        GapSpikeIndices = find((BoutSpikeTimes >= (GapOnset - GapPreMotorLag)) & (BoutSpikeTimes < (GapOffset - GapPreMotorLag)));
        GapFR(INIndex) = length(GapSpikeIndices)/(GapOffset - GapOnset);
        GapNumSpike(INIndex) = length(GapSpikeIndices);
        
        INDur(INIndex) = INOffset - INOnset;
        INGapDur(INIndex) = GapOffset - INOnset;

        FFSpikeIndices = find((BoutSpikeTimes >= (INOnset + FFEdges(1) - FFWindow)) & (BoutSpikeTimes < (INOnset + FFEdges(end) + FFWindow)));
        TempSpikeTimes = BoutSpikeTimes(FFSpikeIndices) - INOnset;
        TempSpikeTimes = TempSpikeTimes(:);
        FFDataStruct(INIndex,:) = histc(TempSpikeTimes, FFDataStructEdges);

        FFSpikeTimes{INIndex} = TempSpikeTimes;
        for FFWinIndex = 1:length(FFEdges),
           FFPST(INIndex, FFWinIndex) = length(find((FFSpikeTimes{INIndex} >= (FFEdges(FFWinIndex) - FFWindow/2)) & (FFSpikeTimes{INIndex} < (FFEdges(FFWinIndex) + FFWindow/2))));
        end
    end
end

MinINGapDur = min(INGapDur);

MaxINs = max([max(Neural_INR.NoofINs) max(Neural_INR.WithinBoutNoofINs(:,1))]);

CorrFs = 10000;
CorrTime = -GapPreMotorLag:1/CorrFs:(min(INGapDur)-GapPreMotorLag);

for i = 1:MaxINs,
    TrialNos(i) = length(find((NoofINs == i) & (Position(:,1) == -1)));
    if (TrialNos(i) < MinNumber)
        continue;
    end
    
    EdgeIndices = find((Edges >= -GapPreMotorLag) & (Edges <= (min(INGapDur) - GapPreMotorLag)));
    
    for j = 1:i,
        if (j == 1)
            FirstINIndices = find((Position(:,2) == 1) & (NoofINs ~=i));
        else
            FirstINIndices = find((Position(:,2) == 1));
        end
        clear FirstINPST LastINPST SpecificINPST;
        FirstINPST = mean(PST(FirstINIndices, EdgeIndices));
        
        if ((j-i-1) == ComparisonIN)
            LastINIndices = find((Position(:,1) == ComparisonIN) & (NoofINs ~=i));
        else
            LastINIndices = find((Position(:,1) == ComparisonIN));
        end
        if (length(LastINIndices) > 1)
            LastINPST = mean(PST(LastINIndices, EdgeIndices));
        else
            LastINPST = (PST(LastINIndices, EdgeIndices));
        end
        
        SpecificINs = find((NoofINs == i) & (Position(:,1) == (j - i - 1)));
        SpecificINPST = mean(PST(SpecificINs, EdgeIndices));
 
        TempCorr = corrcoef((LastINPST - mean(LastINPST)), (SpecificINPST - mean(SpecificINPST)));
        LastINCorr{i}(j) = TempCorr(1,2);
        
        TempCorr = corrcoef((FirstINPST - mean(FirstINPST)), (SpecificINPST - mean(SpecificINPST)));
        FirstINCorr{i}(j) = TempCorr(1,2);
    end
end

for i = 1:MaxINs,
    if (TrialNos(i) < MinPlotNumber)
        continue;
    end
    
    for j = 1:i,
        INGapRaster{i}{j} = [];
        
        SpecificINs = find((NoofINs == i) & (Position(:,1) == (j - i - 1)));
        TrialGapDurs = INGapDur(SpecificINs);
        [SortedGapDurs, SortedIndices] = sort(TrialGapDurs);
        SpecificINs = SpecificINs(SortedIndices);
        INDurs{i}{j} = INDur(SpecificINs);           
        GapDurs{i}{j} = GapDur(SpecificINs);
        INGapDurs{i}{j} = INGapDur(SpecificINs);
        INGapFRs{i}{j} = INGapFR(SpecificINs);
        INFRs{i}{j} = INFR(SpecificINs);
        GapFRs{i}{j} = GapFR(SpecificINs);
        INGapNumSpikes{i}{j} = INGapNumSpike(SpecificINs);
        INPST{i}{j} = PST(SpecificINs,:);
        
        TempPST = [];
        for k = 1:length(SpecificINs),
           INGapRaster{i}{j} = [INGapRaster{i}{j}; [RasterSpikeTimes{SpecificINs(k)} ones(size(RasterSpikeTimes{SpecificINs(k)}))*k]];
           for FFWinIndex = 1:length(FFEdges),
               TempPST(k, FFWinIndex) = length(find((FFSpikeTimes{SpecificINs(k)} >= (FFEdges(FFWinIndex) - FFWindow/2)) & (FFSpikeTimes{SpecificINs(k)} < (FFEdges(FFWinIndex) + FFWindow/2))));
           end
        end
        INGapFF{i}{j} = var(TempPST)./mean(TempPST);
        StartPoint = find(FFEdges <= -GapPreMotorLag, 1, 'last');
        LastGaps = find(Position(:,1) == -1);
        LastGapDurs = INGapDur(LastGaps) - INDur(LastGaps);
        EndPoint = find(FFEdges <= min(LastGapDurs), 1, 'last');
        INGapFFWindowValue{i}(j) = mean(INGapFF{i}{j}(StartPoint:EndPoint));
    end
end

MaxPlots = find(TrialNos >= MinPlotNumber);
MaxPlots = MaxPlots(end);

ValidTrials = find(TrialNos >= MinPlotNumber);
AllFR = ones(length(ValidTrials), MaxPlots)*NaN;
AllFirstINCorr = ones(length(ValidTrials), MaxPlots)*NaN;
AllLastINCorr = ones(length(ValidTrials), MaxPlots)*NaN;
AllFFWinValue = ones(length(ValidTrials), MaxPlots)*NaN;
AllINFR = ones(length(ValidTrials), MaxPlots)*NaN;
AllGapFR = ones(length(ValidTrials), MaxPlots)*NaN;

AllFirstINsCorrValues = [];
for j = 1:length(ValidTrials),
    for k = 1:ValidTrials(j),
        AllFR(j, (MaxPlots - (ValidTrials(j) - k))) = mean(INGapFRs{ValidTrials(j)}{k});
        AllINFR(j, (MaxPlots - (ValidTrials(j) - k))) = mean(INFRs{ValidTrials(j)}{k});
        AllGapFR(j, (MaxPlots - (ValidTrials(j) - k))) = mean(GapFRs{ValidTrials(j)}{k});
        AllFirstINCorr(j, (MaxPlots - (ValidTrials(j) - k))) = FirstINCorr{ValidTrials(j)}(k);
        AllLastINCorr(j, (MaxPlots - (ValidTrials(j) - k))) = LastINCorr{ValidTrials(j)}(k);
        AllFFWinValue(j, (MaxPlots - (ValidTrials(j) - k))) = INGapFFWindowValue{ValidTrials(j)}(k);
        if (k==1)
            AllFirstINsCorrValues = [AllFirstINsCorrValues; FirstINCorr{ValidTrials(j)}(k)];
        end
    end
end

for j = 1:size(AllFR,2),
    TempAllFR(j,:) = [mean(AllFR(find(~isnan(AllFR(:,j))),j)) std(AllFR(find(~isnan(AllFR(:,j))),j))/sqrt(length(find(~isnan(AllFR(:,j))))) j];
    TempAllFirstINCorr(j,:) = [mean(AllFirstINCorr(find(~isnan(AllFirstINCorr(:,j))),j)) std(AllFirstINCorr(find(~isnan(AllFirstINCorr(:,j))),j))/sqrt(length(find(~isnan(AllFirstINCorr(:,j)))))];
    TempAllLastINCorr(j,:) = [mean(AllLastINCorr(find(~isnan(AllLastINCorr(:,j))),j)) std(AllLastINCorr(find(~isnan(AllLastINCorr(:,j))),j))/sqrt(length(find(~isnan(AllLastINCorr(:,j)))))];
    TempAllFFWinValue(j,:) = [mean(AllFFWinValue(find(~isnan(AllFFWinValue(:,j))),j)) std(AllFFWinValue(find(~isnan(AllFFWinValue(:,j))),j))/sqrt(length(find(~isnan(AllFFWinValue(:,j)))))];
    TempAllINFR(j,:) = [mean(AllINFR(find(~isnan(AllINFR(:,j))),j)) std(AllINFR(find(~isnan(AllINFR(:,j))),j))/sqrt(length(find(~isnan(AllINFR(:,j))))) j];
    TempAllGapFR(j,:) = [mean(AllGapFR(find(~isnan(AllGapFR(:,j))),j)) std(AllGapFR(find(~isnan(AllGapFR(:,j))),j))/sqrt(length(find(~isnan(AllGapFR(:,j))))) j];
end        

if (strfind(PlotOption, 'on'))
    
    RasterFig = figure;
    set(gcf, 'Color', 'w');

    FRFig = figure;
    set(gcf, 'Color', 'w');

    CorrFig = figure;
    set(gcf, 'Color', 'w');

    FFFig1 = figure;
    set(gcf, 'Color', 'w');

    AllFRFig = figure;
    set(gcf, 'Color', 'w');

    AllCorrFig = figure;
    set(gcf, 'Color', 'w');

    Rows = length(find(TrialNos >= MinPlotNumber));

    Colours = ['rbmkcg'];
    Index = MaxPlots;
    RowIndex = Rows + 1;


    Index = 0;
    for i = 1:MaxINs,
        PlotHGap = 0.03;
        PlotWGap = 0.035;

        PlotHt = 0.92/Rows - PlotHGap;

        PlotWidth = 0.96/MaxPlots - PlotWGap;

        if (TrialNos(i) >= MinPlotNumber)
            Index = Index + 1;
            RowIndex = RowIndex - 1;
            TempFR = zeros(i,2);
            TempNumSpikes = zeros(i,2);
            for j = 1:i,
                figure(RasterFig);
                TempFR(j,:) = [mean(INGapFRs{i}{j}) std(INGapFRs{i}{j})];
                TempNumSpikes(j,:) = [mean(INGapNumSpikes{i}{j}) std(INGapNumSpikes{i}{j})];

                subplot('Position', [(1 - (i - j + 1)*(PlotWidth + PlotWGap)) (1 - (Index)*(PlotHt + PlotHGap)) PlotWidth PlotHt]);
                PlotRaster(INGapRaster{i}{j}, 'k', 0.25, 0, MinPlotNumber);
                hold on;
%                fill([(ones(1, MinPlotNumber)*0 - GapPreMotorLag) fliplr(INDurs{i}{j}(1:MinPlotNumber) - GapPreMotorLag)], [linspace(0.5, MinPlotNumber+0.5, MinPlotNumber) fliplr(linspace(0.5, MinPlotNumber+0.5, MinPlotNumber))], 'k', 'FaceColor', 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                fill([(ones(1, MinPlotNumber)*0) fliplr(INDurs{i}{j}(1:MinPlotNumber))], [linspace(0.5, MinPlotNumber+0.5, MinPlotNumber) fliplr(linspace(0.5, MinPlotNumber+0.5, MinPlotNumber))], 'k', 'FaceColor', 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                if (j ~= i)
                    fill([(0 + INGapDurs{i}{j}(1:MinPlotNumber)) fliplr(INGapDurs{i}{j}(1:MinPlotNumber) + INDurs{i}{j+1}(1:MinPlotNumber))], [linspace(0.5, MinPlotNumber+0.5, MinPlotNumber) fliplr(linspace(0.5, MinPlotNumber+0.5, MinPlotNumber))], 'k', 'FaceColor', 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                else
                    fill([(0 + INGapDurs{i}{j}(1:MinPlotNumber)) fliplr(INGapDurs{i}{j}(1:MinPlotNumber) + SyllDur(1:MinPlotNumber)')], [linspace(0.5, MinPlotNumber+0.5, MinPlotNumber) fliplr(linspace(0.5, MinPlotNumber+0.5, MinPlotNumber))], 'k', 'FaceColor', 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                end
                if (j ~= 1)
                    fill([(0 - INGapDurs{i}{j-1}(1:MinPlotNumber)) fliplr(0 - INGapDurs{i}{j-1}(1:MinPlotNumber) + INDurs{i}{j-1}(1:MinPlotNumber))], [linspace(0.5, MinPlotNumber+0.5, MinPlotNumber) fliplr(linspace(0.5, MinPlotNumber+0.5, MinPlotNumber))], 'k', 'FaceColor', 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                end
%               fill([(INDurs{i}{j}(1:MinPlotNumber) - GapPreMotorLag) fliplr(INGapDurs{i}{j}(1:MinPlotNumber) - GapPreMotorLag)], [linspace(0.5, MinPlotNumber+0.5, MinPlotNumber) fliplr(linspace(0.5, MinPlotNumber+0.5, MinPlotNumber))], 'k', 'FaceColor', 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                %fill([(INDurs{i}{j}(1:MinPlotNumber)) fliplr(INGapDurs{i}{j}(1:MinPlotNumber))], [linspace(0.5, MinPlotNumber+0.5, MinPlotNumber) fliplr(linspace(0.5, MinPlotNumber+0.5, MinPlotNumber))], 'k', 'FaceColor', 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                axis([-PreTime PostTime -0.5 MinPlotNumber+0.5]);
                plot(Edges(1:end-1)+BinSize/2, mean(INPST{i}{j}(:,1:end-1))/max(max(PST)) * MinPlotNumber, 'Color', 'r', 'LineWidth', 1.5);
                %plot(FFEdges, INGapFF{i}{j}/3 * MinPlotNumber, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
                if ((i == (MaxPlots)) && (j == 1))
                    plot([-PreTime -PreTime], [(MinPlotNumber - 100/max(max(PST))*MinPlotNumber) MinPlotNumber], 'Color', 'r', 'LineWidth', 2);
                end
                %if ((i == (MaxPlots - 1)) && (j == 1))
                %    plot([-PreTime -PreTime], [MinPlotNumber (MinPlotNumber - 1/3*MinPlotNumber)], 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
                %end

                %plot([0 0], [-0.5 MinPlotNumber+0.5], 'k--');
                set(gca, 'Box', 'off');
                set(gca, 'YColor', 'w');
                set(gca, 'TickDir', 'out');
                set(gca, 'XTick', [-PreTime:0.1:PostTime]);
                set(gca, 'TickLen', [0.05 0.025]);

                if (Index ~= MaxPlots)
                    set(gca, 'XTickLabel', []);
                end

                if (i == MaxPlots)
                    figure(FFFig1);
                    set(gcf, 'Color', 'w');
                    subplot('Position', [(1 - (i - j + 1)*(PlotWidth + PlotWGap)) 0.2 PlotWidth 0.75]);
                    Indices = find(Position(:,1) == (j - MaxPlots - 1));
                    [ax{j}, h1, h2] = plotyy(Edges(1:end-1) + BinSize/2, mean(PST(Indices,1:end-1)), FFEdges, var(FFPST(Indices,:))./mean(FFPST(Indices,:)));
                    hold on;
                    axes(ax{j}(1));
                    fill([(Edges(1:end-1)+BinSize/2) fliplr(Edges(1:end-1)+BinSize/2)], [(mean(PST(Indices,1:end-1)) - std(PST(Indices,1:end-1))) fliplr(mean(PST(Indices,1:end-1)) + std(PST(Indices,1:end-1)))], 'k', 'EdgeColor', 'none', 'FaceColor', [1 0.5 0.5], 'FaceAlpha', 0.3);
                    set(gca, 'YColor', 'w');
                    set(gca, 'Box', 'off');
                    set(gca, 'TickDir', 'out');
                    axis tight;
                    tempFRaxis(j,:) = axis;

                    axes(ax{j}(2));
                    set(gca, 'YColor', 'w');
                    set(gca, 'Box', 'off');
                    set(gca, 'TickDir', 'out');
                    set(gca, 'TickLen', [0.03 0.03]);
                    axis tight;
                    tempFFaxis(j,:) = axis;

                    set(h1, 'Color', 'r')
                    set(h2, 'Color', [0.5 0.5 0.5]);
                end
            end

            if (i == MaxPlots)
                for j = 1:i,
                    axes(ax{j}(1));
                    axis([FFEdges(1) FFEdges(end) min(tempFRaxis(:,3)) max(tempFRaxis(:,4))]); 
                    plot([0 0], [min(tempFRaxis(:,3)) max(tempFRaxis(:,4))], 'k--');
                    fill([-GapPreMotorLag min(INGapDur) min(INGapDur) -GapPreMotorLag], [min(tempFRaxis(:,3)) min(tempFRaxis(:,3)) max(tempFRaxis(:,4)) max(tempFRaxis(:,4))], 'k', 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.3);
                    set(gca, 'XTickLabel', []);

                    axes(ax{j}(2));
                    axis([FFEdges(1) FFEdges(end) min(tempFFaxis(:,3)) max(tempFFaxis(:,4))]);
                    set(gca, 'XTickLabel', []);
                end
            end

            figure(FRFig);
            FRFigAxis(Index) = subplot('Position', [0.1 (1 - (Index)*(PlotHt + PlotHGap)) 0.8 PlotHt]);
            %[Ax{Index}, H1, H2] = plotyy([-i:1:-1], TempFR(:,1), [-i:1:-1], INGapFFWindowValue{i});
%            plot([-i:1:-1], TempFR(:,1), 'rs-');
            %axes(Ax{Index}(1));
            hold on;
            set(gca, 'YColor', 'r');
%            set(H1, 'Color', 'r');
            errorbar([-i:1:-1], TempFR(:,1), TempFR(:,2), 'rs-');
            axis([-(MaxPlots+0.25) -0.75 0 max(TempFR(:,1))*1.25]);
            set(gca, 'Box', 'off');
            set(gca, 'YTickLabel', []);
            set(gca, 'XTickLabel', []);
            set(gca, 'TickLen', [0.035 0.025]);

%             axes(Ax{Index}(2));
%             hold on;
%             set(gca, 'YColor', [0.5 0.5 0.5]);
%             set(H2, 'Color', [0.5 0.5 0.5], 'Marker', 'o');
%             axis tight;
%             temp = axis;
%             axis([-(MaxPlots+0.25) -0.75 0 temp(4)]);
%             set(gca, 'Box', 'off');
%             set(gca, 'YTickLabel', []);
%             set(gca, 'XTickLabel', []);
%             set(gca, 'TickLen', [0.035 0.025]);

            figure(CorrFig);
            CorrFigAxis(Index) = subplot('Position', [0.2 (1 - (Index)*(PlotHt + PlotHGap)) 0.75 PlotHt]);
            plot([-i:1:-1], LastINCorr{i}, 'ks-');
            hold on;
%            plot([-i:1:-1], FirstINCorr{i}, 'ks-');
            plot([-(MaxPlots+0.25) -0.75], [0 0], 'k--');
            axis([-(MaxPlots+0.25) -0.75 -1 1]);
            set(gca, 'YTick', [-1 0 1]);
            set(gca, 'TickLen', [0.035 0.025]);
            set(gca, 'Box', 'off');
            set(gca, 'XTickLabel', []);
            set(gca, 'YTickLabel', []);
            
%             if (i == MaxPlots)
%                 figure(AllFRFig);
%                 subplot('Position', [0.1 0.2 0.8 0.75]);
%                 [AllAx, H1, H2] = plotyy(-MaxPlots:1:-1, TempAllFR(:,1), -MaxPlots:1:-1, TempAllFFWinValue(:,1));
% 
%                 axes(AllAx(1));
%                 hold on;
%                 set(gca, 'YColor', 'r');
%                 errorbar(-MaxPlots:1:-1, TempAllFR(:,1), TempAllFR(:,2), 'ro-');
%                 set(gca, 'Box', 'off');
%                 axis tight;
%                 temp = axis;
%                 axis([(-MaxPlots - 0.25) -0.75 0 temp(4)]);
%                 set(gca, 'TickLen', [0.035 0.025]);
%                 set(gca, 'Box', 'off');
%                 set(gca, 'XTickLabel', []);
%                 set(gca, 'YTickLabel', []);
% 
%                 axes(AllAx(2));
%                 hold on;
%                 errorbar(-MaxPlots:1:-1, TempAllFFWinValue(:,1), TempAllFFWinValue(:,2), 'ko-', 'Color', [0.5 0.5 0.5]);
%                 set(gca, 'Box', 'off');
%                 set(gca, 'YColor', [0.5 0.5 0.5]);
%                 axis tight;
%                 temp = axis;
%                 axis([(-MaxPlots - 0.25) -0.75 0 temp(4)]);
%                 axis
%                 set(gca, 'TickLen', [0.035 0.025]);
%                 set(gca, 'Box', 'off');
%                 set(gca, 'XTickLabel', []);
%                 set(gca, 'YTickLabel', []);
% 
%                 figure(AllCorrFig);
%                 errorbar(-MaxPlots:1:-1, TempAllFirstINCorr(:,1), TempAllFirstINCorr(:,2), 'ks-');
%                 hold on;
%                 errorbar(-MaxPlots:1:-1, TempAllLastINCorr(:,1), TempAllLastINCorr(:,2), 'ks-', 'MarkerFaceColor', 'k');
%                 plot([-(MaxPlots+0.25) -0.75], [0 0], 'k--');
% 
%                 set(gca, 'Position', [0.2 0.2 0.75 0.75]);
%                 set(gca, 'Box', 'off');
%                 axis tight;
%                 temp = axis;
%                 axis([(-MaxPlots - 0.25) -0.75 -1 1]);
%                 set(gca, 'YTick', [-1 0 1]);
%                 set(gca, 'TickLen', [0.035 0.025]);
%                 set(gca, 'Box', 'off');
%                 set(gca, 'XTickLabel', []);
%                 set(gca, 'YTickLabel', []);
%             end
        end
    end

    for i = 1:length(FRFigAxis),
%         axes(Ax{i}(1));
        axes(FRFigAxis(i));
        axis tight;
        TempIndFRaxis(i,:) = axis;

%         axes(Ax{i}(2));
%         axis tight;
%         TempIndFFaxis(i,:) = axis;
    end

    TempIndFR(1:3) = [-(MaxPlots + 0.25) -0.75 0];
    TempIndFR(4) = 1.05*max(TempIndFRaxis(:,4));

%     TempIndFF(1:3) = [-(MaxPlots + 0.25) -0.75 0];
%     TempIndFF(4) = 1.05*max(TempIndFFaxis(:,4));

    for i = 1:length(FRFigAxis),
        axes(FRFigAxis(i));
        axis(TempIndFR);
%         axes(Ax{i}(2));
%         axis(TempIndFF);    
    end
end

TempAllFR(:,1) = TempAllFR(:,1) - mean(INGapFR);

LastvsOtherIndex = (mean(INGapFR(find(Position(:,1) == -1))) - mean(INGapFR(find(Position(:,1) < -1))))/(mean(INGapFR(find(Position(:,1) < -1))) + mean(INGapFR(find(Position(:,1) ~= -1))));

LastvsOtherFR(1) = mean(INGapFR(find(Position(:,1) == -1)));
LastvsOtherFR(2) = std(INGapFR(find(Position(:,1) == -1)));
LastvsOtherFR(3) = mean(INGapFR(find(Position(:,1) < -1)));
LastvsOtherFR(4) = std(INGapFR(find(Position(:,1) < -1)));

Indices = find(Position(:,1) <= -1);
Indices = Indices(:);
LastvsOtherPosFR = [Position(Indices,1) INGapFR(Indices)'];
for i = min(LastvsOtherPosFR(:,1)):1:max(LastvsOtherPosFR(:,1)),
    if (length(find(LastvsOtherPosFR(:,1) == i) >= 2))
        break;
    end
end
Indices = find(LastvsOtherPosFR(:,1) >=i);
[LastvsOthersSig, anovatab, stats] = kruskalwallis(LastvsOtherPosFR(Indices,2), LastvsOtherPosFR(Indices,1), 'off');
figure;
multcompare(stats);

LastvsOtherPosFR(find(LastvsOtherPosFR(:,1) < -1), 1) = -2;
AllSig = kruskalwallis(LastvsOtherPosFR(:,2), LastvsOtherPosFR(:,1), 'off');

% LastINFR = INGapFR(find(Position == -1));
% OtherINFR = INGapFR(find(Position < -1));
% Temp = [];
% for j = 1:length(OtherINFR),
%     Temp = [Temp; (LastINFR - OtherINFR(j))];
% end
% Temp = abs(Temp(:));
% LastvsOtherFRDiff(1,:) = [mean(Temp) std(Temp)/sqrt(length(LastINFR))];
% Temp2 = [];
% for j = 1:length(OtherINFR),
%     for k = j+1:length(OtherINFR),
%         Temp2 = [Temp2; OtherINFR(j) - OtherINFR(k)];
%     end
% end
% Temp2 = abs(Temp2(:));
% OthervsOtherFRDiff(1,:) = [mean(Temp2) std(Temp2)/sqrt(length(LastINFR))];
% Sig = anova1([Temp; Temp2], [ones(size(Temp))*1; ones(size(Temp2))*2]);    


disp('Finished');