function [INNeuralAnalysisResults] = IntroNoteNeuralAnalysisINSequence(Neural_INR, PlotOption, MinTrialNumber, PlotPSTH, varargin)

% Using a system for AllINPosition that keeps track of whether the intro
% note was the first, last or a middle intro note. The way I do this is by
% having 3 boolean flags for first, middle and last intro note. For
% instance if there was only one intro note, it would have the flags 1 0 1
% to indicate that is the first, it is also the last. 

% It also has 3 more flags - one for the total number of INs in that
% sequence, the second one gives the trial # for the particular sequence of
% INs and the third gives the position of the IN within that sequence. With
% all of this it should be easy to reconstruct the position and trial no of
% each IN.

Width = 0.01;
GaussianLen = 4;
IFRFs = 2000;
XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (IFRFs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * IFRFs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * IFRFs) * (Width * IFRFs)));


AllINSpikeTrain = [];
AllINPosition = [];
AllINIndex = 0;
AllINIFR = [];
AllINFR = [];
AllINPST = [];

MinNumber = 2;
if (nargin > 4),
    BinSize = varargin{1};
else
    BinSize = 0.0025;
end

if (nargin > 5)
    Edges = varargin{2}:BinSize:varargin{3};
else
    Edges = -0.1:BinSize:0.1;
end

if (nargin > 7)
    PlotMinINs = varargin{4};
end

IFREdges = Edges(1):1/IFRFs:Edges(end);

MaxINs = max([max(Neural_INR.NoofINs) max(Neural_INR.WithinBoutNoofINs(:,1))]);

for i = 1:MaxINs,
    for j = 1:i+1,
        INIFR{i}{j} = [];
        INRaster{i}{j} = [];
        PrevSyllOffsetRaster{i}{j} = [];
        SyllOffsetRaster{i}{j} = [];
        NextSyllOnsetRaster{i}{j} = [];
    end
end

for i = 1:MaxINs,
    TrialNo = 0;
    for j = 1:length(Neural_INR.NoofINs),
        if (Neural_INR.NoofINs(j) == i)
            TrialNo = TrialNo + 1;
            StartTimes = [Neural_INR.BoutDetails(j).onsets(Neural_INR.INs{j}); Neural_INR.BoutDetails(j).onsets(Neural_INR.INs{j}(end)+1)];
            EndTimes = [Neural_INR.BoutDetails(j).offsets(Neural_INR.INs{j}); Neural_INR.BoutDetails(j).offsets(Neural_INR.INs{j}(end)+1)];
            
            if (Neural_INR.INs{j}(1) ~= 1)
                PrevSyllOffsetTimes = [Neural_INR.BoutDetails(j).offsets(Neural_INR.INs{j}-1); Neural_INR.BoutDetails(j).offsets(Neural_INR.INs{j}(end))];
            else
                if (i > 1)
                    PrevSyllOffsetTimes = [NaN; Neural_INR.BoutDetails(j).offsets(Neural_INR.INs{j}(2:end) - 1); Neural_INR.BoutDetails(j).offsets(Neural_INR.INs{j}(end))];
                else
                    PrevSyllOffsetTimes = [NaN; Neural_INR.BoutDetails(j).offsets(Neural_INR.INs{j}(end))];
                end
            end
            SyllOffsetTimes = [Neural_INR.BoutDetails(j).offsets(Neural_INR.INs{j}); Neural_INR.BoutDetails(j).offsets(Neural_INR.INs{j}(end) + 1)];
            
            if ((Neural_INR.INs{j}(end) + 2) > length(Neural_INR.BoutDetails(j).onsets))
                NextSyllOnsetTimes = [Neural_INR.BoutDetails(j).onsets(Neural_INR.INs{j}+1); NaN];
            else
                NextSyllOnsetTimes = [Neural_INR.BoutDetails(j).onsets(Neural_INR.INs{j}+1); Neural_INR.BoutDetails(j).offsets(Neural_INR.INs{j}(end) + 2)];
            end
            
            for k = 1:i+1,
                IFRStartIndex = find(Neural_INR.BoutDetails(j).IFR(1,:) <= (StartTimes(k) + Edges(1)), 1, 'last');
                IFREndIndex = find(Neural_INR.BoutDetails(j).IFR(1,:) <= (StartTimes(k) + Edges(end)), 1, 'last');
                INIFR{i}{k} = [INIFR{i}{k}; interp1(Neural_INR.BoutDetails(j).IFR(1, IFRStartIndex:IFREndIndex), Neural_INR.BoutDetails(j).IFR(2, IFRStartIndex:IFREndIndex), StartTimes(k) + IFREdges)];
                
                PreMotorSpikeIndices = find((Neural_INR.BoutDetails(j).SpikeTimes >= (StartTimes(k)-0.045)) & (Neural_INR.BoutDetails(j).SpikeTimes < (EndTimes(k) - 0.045)));
                
                SpikeIndices = find((Neural_INR.BoutDetails(j).SpikeTimes >= (StartTimes(k)+Edges(1))) & (Neural_INR.BoutDetails(j).SpikeTimes < (StartTimes(k)+Edges(end))));
                if (~isempty(SpikeIndices))
                    INRaster{i}{k} = [INRaster{i}{k}; [(Neural_INR.BoutDetails(j).SpikeTimes(SpikeIndices) - StartTimes(k)) ones(length(SpikeIndices),1)*TrialNo]];
                end
                INFR{i}(TrialNo, k) = length(SpikeIndices)/(Edges(end) - Edges(1));
                INFR_IFR{i}(TrialNo, k) = mean(INIFR{i}{k}(end,1:end-1));
                INPST{i}{k}(TrialNo,:) = histc(Neural_INR.BoutDetails(j).SpikeTimes, StartTimes(k) + Edges)/BinSize;
                PrevSyllOffsetRaster{i}{k} = [PrevSyllOffsetRaster{i}{k}; [(PrevSyllOffsetTimes(k)-StartTimes(k)) TrialNo]];
                SyllOffsetRaster{i}{k} = [SyllOffsetRaster{i}{k}; [(SyllOffsetTimes(k)-StartTimes(k)) TrialNo]];
                NextSyllOnsetRaster{i}{k} = [NextSyllOnsetRaster{i}{k}; [(NextSyllOnsetTimes(k)-StartTimes(k)) TrialNo]];
                if (k < (i+1))
                    AllINIndex = AllINIndex + 1;
                    
                    AllINPreMotorFR(AllINIndex,:) = length(PreMotorSpikeIndices)/(EndTimes(k) - StartTimes(k));
                    AllINPreMotorNumSpikes(AllINIndex,:) = length(PreMotorSpikeIndices);
                    if (~isempty(PreMotorSpikeIndices))
                        AllINPreMotorSpikeTrain{AllINIndex,:} = (Neural_INR.BoutDetails(j).SpikeTimes(PreMotorSpikeIndices) - (StartTimes(k) - 0.045));
                    else
                        AllINPreMotorSpikeTrain{AllINIndex,:} = [];
                    end
                    
                    IntervalPreMotorSpikeIndices = find((Neural_INR.BoutDetails(j).SpikeTimes >= (EndTimes(k)-0.045)) & (Neural_INR.BoutDetails(j).SpikeTimes < (StartTimes(k+1) - 0.045)));
                    AllINIntPreMotorFR(AllINIndex,:) = length(IntervalPreMotorSpikeIndices)/(StartTimes(k+1) - EndTimes(k));
                    AllINIntLen(AllINIndex,:) = StartTimes(k+1) - EndTimes(k);
                    AllINIntPreMotorNumSpikes(AllINIndex,:) = length(IntervalPreMotorSpikeIndices);
                    if (~isempty(IntervalPreMotorSpikeIndices))
                        AllINIntPreMotorSpikeTrain{AllINIndex,:} = (Neural_INR.BoutDetails(j).SpikeTimes(IntervalPreMotorSpikeIndices) - (EndTimes(k) - 0.045));
                    else
                        AllINIntPreMotorSpikeTrain{AllINIndex,:} = [];
                    end
                    
                    AllINFeats(AllINIndex,:) = Neural_INR.BoutDetails(j).Feats(Neural_INR.INs{j}(k),1:4);
                    AllINFR(AllINIndex,:) = length(SpikeIndices)/(Edges(end) - Edges(1));
                    AllINIFR(AllINIndex,:) = interp1(Neural_INR.BoutDetails(j).IFR(1, IFRStartIndex:IFREndIndex), Neural_INR.BoutDetails(j).IFR(2, IFRStartIndex:IFREndIndex), StartTimes(k) + IFREdges);
                    AllINFR_IFR(AllINIndex,:) = mean(AllINIFR(AllINIndex,1:end-1));
                    if (~isempty(SpikeIndices))
                        AllINSpikeTrain{AllINIndex} = (Neural_INR.BoutDetails(j).SpikeTimes(SpikeIndices) - StartTimes(k));
                        AllINPST(AllINIndex,:) = histc(AllINSpikeTrain{AllINIndex}, Edges)/BinSize;
                    else
                        AllINSpikeTrain{AllINIndex} = [];
                        AllINPST(AllINIndex,:) = zeros(size(Edges));
                    end
                    if ((k == 1) && (i == 1))
                        AllINPosition(AllINIndex,:) = [1 0 1 i TrialNo k (k-i-1)];
                    else
                        if (k == 1)
                            AllINPosition(AllINIndex,:) = [1 0 0 i TrialNo k (k-i-1)];
                        else
                            if (k == i)
                                AllINPosition(AllINIndex,:) = [0 0 1 i TrialNo k (k-i-1)];
                            else
                                AllINPosition(AllINIndex,:) = [0 1 0 i TrialNo k (k-i-1)];
                            end
                        end
                    end
                end
            end
        end
    end
end

      
for i = 1:MaxINs,
    if (isempty(SyllOffsetRaster{i}{1}))
        TrialNo = 0;
    else
        TrialNo = SyllOffsetRaster{i}{1}(end,2);
    end
    for j = 1:size(Neural_INR.WithinBoutNoofINs,1),
        if (Neural_INR.WithinBoutNoofINs(j,1) == i)
            TrialNo = TrialNo + 1;
            StartTimes = [Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).onsets(Neural_INR.WithinBoutINs{j}); Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).onsets(Neural_INR.WithinBoutINs{j}(end)+1)];
            EndTimes = [Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).offsets(Neural_INR.WithinBoutINs{j}); Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).offsets(Neural_INR.WithinBoutINs{j}(end)+1)];
            
            if (Neural_INR.WithinBoutINs{j}(1) ~= 1)
                PrevSyllOffsetTimes = [Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).offsets(Neural_INR.WithinBoutINs{j}-1); Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).offsets(Neural_INR.WithinBoutINs{j}(end))];
            else
                if (i > 1)
                    PrevSyllOffsetTimes = [NaN; Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).offsets(Neural_INR.WithinBoutINs{j}(2:end) - 1); Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).offsets(Neural_INR.WithinBoutINs{j}(end))];
                else
                    PrevSyllOffsetTimes = [NaN; Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).offsets(Neural_INR.WithinBoutINs{j}(end))];
                end
            end
            SyllOffsetTimes = [Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).offsets(Neural_INR.WithinBoutINs{j}); Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).offsets(Neural_INR.WithinBoutINs{j}(end) + 1)];
            
            if ((Neural_INR.WithinBoutINs{j}(end) + 2) > length(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).onsets))
                NextSyllOnsetTimes = [Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).onsets(Neural_INR.WithinBoutINs{j}+1); NaN];
            else
                NextSyllOnsetTimes = [Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).onsets(Neural_INR.WithinBoutINs{j}+1); Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).offsets(Neural_INR.WithinBoutINs{j}(end) + 2)];
            end
            
            for k = 1:i+1,
                
                IFRStartIndex = find(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).IFR(1,:) <= (StartTimes(k) + Edges(1)), 1, 'last');
                IFREndIndex = find(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).IFR(1,:) <= (StartTimes(k) + Edges(end)), 1, 'last');
                INIFR{i}{k} = [INIFR{i}{k}; interp1(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).IFR(1, IFRStartIndex:IFREndIndex), Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).IFR(2, IFRStartIndex:IFREndIndex), StartTimes(k) + IFREdges)];
                
                PreMotorSpikeIndices = find((Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).SpikeTimes >= (StartTimes(k) - 0.045)) & (Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).SpikeTimes < (EndTimes(k) - 0.045)));
                
                SpikeIndices = find((Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).SpikeTimes >= (StartTimes(k)+Edges(1))) & (Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).SpikeTimes < (StartTimes(k)+Edges(end))));
                if (~isempty(SpikeIndices))
                    INRaster{i}{k} = [INRaster{i}{k}; [(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).SpikeTimes(SpikeIndices) - StartTimes(k)) ones(length(SpikeIndices),1)*TrialNo]];
                end
                INFR{i}(TrialNo, k) = length(SpikeIndices)/(Edges(end) - Edges(1));
                INFR_IFR{i}(TrialNo, k) = mean(INIFR{i}{k}(end,1:end-1));
                INPST{i}{k}(TrialNo,:) = histc(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).SpikeTimes, StartTimes(k) + Edges)/BinSize;
                PrevSyllOffsetRaster{i}{k} = [PrevSyllOffsetRaster{i}{k}; [(PrevSyllOffsetTimes(k)-StartTimes(k)) TrialNo]];
                SyllOffsetRaster{i}{k} = [SyllOffsetRaster{i}{k}; [(SyllOffsetTimes(k)-StartTimes(k)) TrialNo]];
                NextSyllOnsetRaster{i}{k} = [NextSyllOnsetRaster{i}{k}; [(NextSyllOnsetTimes(k)-StartTimes(k)) TrialNo]];
                if (k < (i+1))
                    AllINIndex = AllINIndex + 1;
                    AllINPreMotorFR(AllINIndex,:) = length(PreMotorSpikeIndices)/(EndTimes(k) - StartTimes(k));
                    AllINPreMotorNumSpikes(AllINIndex,:) = length(PreMotorSpikeIndices);
                    if (~isempty(PreMotorSpikeIndices))
                        AllINPreMotorSpikeTrain{AllINIndex,:} = (Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).SpikeTimes(SpikeIndices) - (StartTimes(k) - 0.045));
                    else
                        AllINPreMotorSpikeTrain{AllINIndex,:} = [];
                    end
                    
                    IntervalPreMotorSpikeIndices = find((Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).SpikeTimes >= (EndTimes(k) - 0.045)) & (Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).SpikeTimes < (StartTimes(k+1) - 0.045)));
                    AllINIntPreMotorFR(AllINIndex,:) = length(IntervalPreMotorSpikeIndices)/(StartTimes(k+1) - EndTimes(k));
                    AllINIntLen(AllINIndex,:) = StartTimes(k+1) - EndTimes(k);
                    AllINIntPreMotorNumSpikes(AllINIndex,:) = length(IntervalPreMotorSpikeIndices);
                    
                    if (~isempty(IntervalPreMotorSpikeIndices))
                        AllINIntPreMotorSpikeTrain{AllINIndex,:} = (Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).SpikeTimes(SpikeIndices) - (EndTimes(k) - 0.045));
                    else
                        AllINIntPreMotorSpikeTrain{AllINIndex,:} = [];
                    end
                    
                    AllINFeats(AllINIndex,:) = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).Feats(Neural_INR.WithinBoutINs{j}(k),1:4);
                    AllINFR(AllINIndex,:) = length(SpikeIndices)/(Edges(end) - Edges(1));
                    AllINIFR(AllINIndex,:) = interp1(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).IFR(1, IFRStartIndex:IFREndIndex), Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).IFR(2, IFRStartIndex:IFREndIndex), StartTimes(k) + IFREdges);
                    AllINFR_IFR(AllINIndex,:) = mean(AllINIFR(AllINIndex,1:end-1));
                    if (~isempty(SpikeIndices))
                        AllINSpikeTrain{AllINIndex} = (Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).SpikeTimes(SpikeIndices) - StartTimes(k));
                        AllINPST(AllINIndex,:) = histc(AllINSpikeTrain{AllINIndex}, Edges)/BinSize;
                    else
                        AllINSpikeTrain{AllINIndex} = [];
                        AllINPST(AllINIndex,:) = zeros(size(Edges));
                    end
                    if ((k == 1) && (i == 1))
                        AllINPosition(AllINIndex,:) = [1 0 1  i TrialNo k (k-i-1)];
                    else
                        if (k == 1)
                            AllINPosition(AllINIndex,:) = [1 0 0  i TrialNo k (k-i-1)];
                        else
                            if (k == i)
                                AllINPosition(AllINIndex,:) = [0 0 1  i TrialNo k (k-i-1)];
                            else
                                AllINPosition(AllINIndex,:) = [0 1 0  i TrialNo k (k-i-1)];
                            end
                        end
                    end
                end
            end
        end
    end
end

for i = 1:MaxINs,
    if (isempty(SyllOffsetRaster{i}{1}))
        NumTrials(i) = 0;
    else
        NumTrials(i) = max(SyllOffsetRaster{i}{1}(:,2));
    end
end
disp(NumTrials);

Indices = find(NumTrials > MinTrialNumber);
MinTrials = min(NumTrials(Indices));
%if (MinTrials > 10)
%    MinTrials = 10;
%end

Temp_INIFR = INIFR;

% now do the within position and across position distances
ValidIntroNoteNos = [];
INNeuralAnalysisResults.AllINFR_Ratios = [];
for i = 1:length(INFR),
    if (i > 1)
        if (~isempty(INFR{i}))
            TempRatios = INFR{i}(:,2:end-1)./INFR{i}(:,1:end-2);
            TempRatios = reshape(TempRatios, size(TempRatios,1)*size(TempRatios,2), 1);
            INNeuralAnalysisResults.AllINFR_Ratios = [INNeuralAnalysisResults.AllINFR_Ratios; TempRatios];
        end
    end
    if (size(INFR{i},1) >= 2)
        ValidIntroNoteNos = [ValidIntroNoteNos; i];
    end
end

INNeuralAnalysisResults.ValidIntroNoteNos = ValidIntroNoteNos;

SeqIndex = 1;
INNeuralAnalysisResults.CommonLast.SimilarPositionFRDistances = [];
INNeuralAnalysisResults.CommonLast.OtherPositionFRDistances = [];
INNeuralAnalysisResults.CommonLast.SimilarPositionCorr = [];
INNeuralAnalysisResults.CommonLast.OtherPositionCorr = [];
for i = ValidIntroNoteNos',
    Indices1 = find(AllINPosition(:,4) == i);
    for j = min(AllINPosition(Indices1,end)):1:max(AllINPosition(Indices1,end)),
        TempSamePosDist = [];
        TempOtherPosDist = [];
        TempSamePosCorr = [];
        TempOtherPosCorr = [];
        INGroup1 = find((AllINPosition(:,end) == j) & (AllINPosition(:,4) == i));
        for k = ValidIntroNoteNos',
            Indices2 = find(AllINPosition(:,4) == k);
            for m = min(AllINPosition(Indices2,end)):1:max(AllINPosition(Indices2,end)),
                INGroup2 = find((AllINPosition(:,end) == m) & (AllINPosition(:,4) == k));
                if ((j == m) && (i == k))
                    continue;
                else
                    FRDistance = abs(mean(AllINFR(INGroup1)) - mean(AllINFR(INGroup2)));
                    
                    INGroup1IFR = conv(mean(AllINIFR(INGroup1,1:end-1)), GaussWin, 'same');
                    INGroup1IFR = INGroup1IFR - mean(INGroup1IFR);
                    
                    INGroup2IFR = conv(mean(AllINIFR(INGroup2,1:end-1)), GaussWin, 'same');
                    INGroup2IFR = INGroup2IFR - mean(INGroup2IFR);
                    
                    Corr = (INGroup1IFR*INGroup2IFR')./(norm(INGroup1IFR)*norm(INGroup2IFR));
                    
                    if (j == m)
                        TempSamePosDist = [TempSamePosDist; FRDistance];
                        TempSamePosCorr = [TempSamePosCorr; Corr];
                        INNeuralAnalysisResults.CommonLast.SimilarPositionFRDistances = [INNeuralAnalysisResults.CommonLast.SimilarPositionFRDistances; FRDistance];
                        INNeuralAnalysisResults.CommonLast.SimilarPositionCorr = [INNeuralAnalysisResults.CommonLast.SimilarPositionCorr; Corr];
                    else
                        TempOtherPosDist = [TempOtherPosDist; FRDistance];
                        TempOtherPosCorr = [TempOtherPosCorr; Corr];
                        INNeuralAnalysisResults.CommonLast.OtherPositionFRDistances = [INNeuralAnalysisResults.CommonLast.OtherPositionFRDistances; FRDistance];
                        INNeuralAnalysisResults.CommonLast.OtherPositionCorr = [INNeuralAnalysisResults.CommonLast.OtherPositionCorr; Corr];
                    end
                end
            end
        end
        if (~isempty(TempSamePosDist) & ~isempty(TempOtherPosDist))
            INNeuralAnalysisResults.CommonLastIndividualDistances(SeqIndex, :) = [mean(TempSamePosDist) mean(TempOtherPosDist) i j];
            INNeuralAnalysisResults.CommonLastIndividualCorrs(SeqIndex, :) = [mean(TempSamePosCorr) mean(TempOtherPosCorr) i j];
            SeqIndex = SeqIndex + 1;
        end
    end
end

SeqIndex = 1;
INNeuralAnalysisResults.CommonFirst.SimilarPositionFRDistances = [];
INNeuralAnalysisResults.CommonFirst.OtherPositionFRDistances = [];
INNeuralAnalysisResults.CommonFirst.SimilarPositionCorr = [];
INNeuralAnalysisResults.CommonFirst.OtherPositionCorr = [];
for i = ValidIntroNoteNos',
    Indices1 = find(AllINPosition(:,4) == i);
    for j = min(AllINPosition(Indices1,end-1)):1:max(AllINPosition(Indices1,end-1)),
        TempSamePosDist = [];
        TempOtherPosDist = [];
        TempSamePosCorr = [];
        TempOtherPosCorr = [];
        INGroup1 = find((AllINPosition(:,end-1) == j) & (AllINPosition(:,4) == i));
        for k = ValidIntroNoteNos',
            Indices2 = find(AllINPosition(:,4) == k);
            for m = min(AllINPosition(Indices2,end-1)):1:max(AllINPosition(Indices2,end-1)),
                INGroup2 = find((AllINPosition(:,end-1) == m) & (AllINPosition(:,4) == k));
                if ((j == m) && (i == k))
                    continue;
                else
                    FRDistance = abs(mean(AllINFR(INGroup1)) - mean(AllINFR(INGroup2)));
                    
                    INGroup1IFR = conv(mean(AllINIFR(INGroup1,1:end-1)), GaussWin, 'same');
                    INGroup1IFR = INGroup1IFR - mean(INGroup1IFR);
                    
                    INGroup2IFR = conv(mean(AllINIFR(INGroup2,1:end-1)), GaussWin, 'same');
                    INGroup2IFR = INGroup2IFR - mean(INGroup2IFR);
                    
                    Corr = (INGroup1IFR*INGroup2IFR')./(norm(INGroup1IFR)*norm(INGroup2IFR));
                    
                    if (j == m)
                        TempSamePosDist = [TempSamePosDist; FRDistance];
                        TempSamePosCorr = [TempSamePosCorr; Corr];
                        INNeuralAnalysisResults.CommonFirst.SimilarPositionFRDistances = [INNeuralAnalysisResults.CommonFirst.SimilarPositionFRDistances; FRDistance];
                        INNeuralAnalysisResults.CommonFirst.SimilarPositionCorr = [INNeuralAnalysisResults.CommonFirst.SimilarPositionCorr; Corr];
                    else
                        TempOtherPosDist = [TempOtherPosDist; FRDistance];
                        TempOtherPosCorr = [TempOtherPosCorr; Corr];
                        INNeuralAnalysisResults.CommonFirst.OtherPositionFRDistances = [INNeuralAnalysisResults.CommonFirst.OtherPositionFRDistances; FRDistance];
                        INNeuralAnalysisResults.CommonFirst.OtherPositionCorr = [INNeuralAnalysisResults.CommonFirst.OtherPositionCorr; Corr];
                    end
                end
            end
        end
        if (~isempty(TempSamePosDist) & ~isempty(TempOtherPosDist))
            INNeuralAnalysisResults.CommonFirstIndividualDistances(SeqIndex, :) = [mean(TempSamePosDist) mean(TempOtherPosDist) i j];
            INNeuralAnalysisResults.CommonFirstIndividualCorrs(SeqIndex, :) = [mean(TempSamePosCorr) mean(TempOtherPosCorr) i j];
            SeqIndex = SeqIndex + 1;
        end
    end
end

clear Corr;
LastINIndices = find(AllINPosition(:,end) == -1);

for i = 1:MaxINs,
    INIndices = find(AllINPosition(:,4) == i);
    if (length(INIndices)/i > MinTrialNumber)
        for j = 1:1:i,
            INIndices = find((AllINPosition(:,4) == i) & (AllINPosition(:,6) == j));
%            INDistance{i}(:,j) = pdist2((AllINFeats(INIndices,:)), mean(AllINFeats(LastINIndices,:)), 'mahalanobis', cov(AllINFeats(LastINIndices,:)));
        end
    end
end
        
for i = 1:MaxINs,
    for j = 1:length(Temp_INIFR{i}),
        if (size(Temp_INIFR{i}{j},1) > MinTrialNumber)
            Corr{i}(j) = 0;
            TempCorr = [];
            
            for k = 1:size(Temp_INIFR{i}{j},1),
                Temp_INIFR{i}{j}(k,1:end-1) = conv(Temp_INIFR{i}{j}(k,1:end-1), GaussWin, 'same');
            end
            Temp_INIFR{i}{j}(:,1:end-1) = Temp_INIFR{i}{j}(:,1:end-1) - repmat(mean(Temp_INIFR{i}{j}(:,1:end-1), 2), 1, size(Temp_INIFR{i}{j}(:,1:end-1),2));
            for FirstTrain = 1:size(Temp_INIFR{i}{j},1),
                for SecondTrain = FirstTrain+1:size(Temp_INIFR{i}{j},1),
                    TempCorr = [TempCorr; (Temp_INIFR{i}{j}(FirstTrain,1:end-1)*Temp_INIFR{i}{j}(SecondTrain,1:end-1)'/(norm(Temp_INIFR{i}{j}(FirstTrain,1:end-1))*norm(Temp_INIFR{i}{j}(SecondTrain,1:end-1))))];
                end
            end
            Corr{i}(j) = mean(TempCorr);
            StdCorr{i}(j) = std(TempCorr);
        end
    end
end

for i = 1:MaxINs,
    for j = 1:length(INRaster{i}),
        if (~isempty(SyllOffsetRaster{i}{j}))
            [SortedDurations, SortedIndices] = sort(SyllOffsetRaster{i}{j}(:,1));
            for k = 1:length(SortedIndices),
                if (~isempty(INRaster{i}{j}))
                    TempIndices = find(INRaster{i}{j}(:,2) == SortedIndices(k));
                    INRaster{i}{j}(TempIndices,2) = 1000+k;
                end
                TempIndices = find(PrevSyllOffsetRaster{i}{j}(:,2) == SortedIndices(k));
                PrevSyllOffsetRaster{i}{j}(TempIndices,2) = 1000+k;

                TempIndices = find(SyllOffsetRaster{i}{j}(:,2) == SortedIndices(k));
                SyllOffsetRaster{i}{j}(TempIndices,2) = 1000+k;

                TempIndices = find(NextSyllOnsetRaster{i}{j}(:,2) == SortedIndices(k));
                NextSyllOnsetRaster{i}{j}(TempIndices,2) = 1000+k;
            end

            if (~isempty(INRaster{i}{j}))
                INRaster{i}{j}(:,2) = INRaster{i}{j}(:,2) - 1000;
            end

            PrevSyllOffsetRaster{i}{j}(:,2) = PrevSyllOffsetRaster{i}{j}(:,2) - 1000;
            [SortedDurations, SortedIndices] = sort(PrevSyllOffsetRaster{i}{j}(:,2));
            PrevSyllOffsetRaster{i}{j} = PrevSyllOffsetRaster{i}{j}(SortedIndices,:);

            SyllOffsetRaster{i}{j}(:,2) = SyllOffsetRaster{i}{j}(:,2) - 1000;
            [SortedDurations, SortedIndices] = sort(SyllOffsetRaster{i}{j}(:,2));
            SyllOffsetRaster{i}{j} = SyllOffsetRaster{i}{j}(SortedIndices,:);

            NextSyllOnsetRaster{i}{j}(:,2) = NextSyllOnsetRaster{i}{j}(:,2) - 1000;
            [SortedDurations, SortedIndices] = sort(NextSyllOnsetRaster{i}{j}(:,2));
            NextSyllOnsetRaster{i}{j} = NextSyllOnsetRaster{i}{j}(SortedIndices,:);
        end
    end
end



FirstINs = find(AllINPosition(:,1) == 1);
FirstINST = AllINIFR(FirstINs,1:end-1);
for i = 1:size(FirstINST, 1),
    FirstINST(i,:) = conv(FirstINST(i,:), GaussWin, 'same');
end
FirstINST = FirstINST - repmat(mean(FirstINST, 2), 1, size(FirstINST,2));

NormFirstINST = [];
for i = 1:size(FirstINST,1),
    NormFirstINST(i) = norm(FirstINST(i,:));
end
TempCorr = [];
for i = 1:size(FirstINST,1),
    TempCorr = [TempCorr (FirstINST(i,:)*FirstINST(i+1:end,:)')./(NormFirstINST(i)*NormFirstINST(i+1:end))];
end
INNeuralAnalysisResults.FirstINs_Corr = [mean(TempCorr(find(~isnan(TempCorr)))) std(TempCorr(find(~isnan(TempCorr))))];

LastINs = find(AllINPosition(:,3) == 1);
LastINST = AllINIFR(LastINs,1:end-1);
for i = 1:size(LastINST, 1),
    LastINST(i,:) = conv(LastINST(i,:), GaussWin, 'same');
end
LastINST = LastINST - repmat(mean(LastINST, 2), 1, size(LastINST,2));
NormLastINST = [];
for i = 1:size(LastINST,1),
    NormLastINST(i) = norm(LastINST(i,:));
end
TempCorr = [];
for i = 1:size(LastINST,1),
    TempCorr = [TempCorr (LastINST(i,:)*LastINST(i+1:end,:)')./(NormLastINST(i)*NormLastINST(i+1:end))];
end
INNeuralAnalysisResults.LastINs_Corr = [mean(TempCorr(find(~isnan(TempCorr)))) std(TempCorr(find(~isnan(TempCorr))))];

for i = 1:MaxINs,
    for j = 1:length(Temp_INIFR{i}),
        if (size(Temp_INIFR{i}{j},1) > MinTrialNumber)
              
            for k = 1:size(Temp_INIFR{i}{j},1),
                Temp_INIFR{i}{j}(k,1:end-1) = conv(Temp_INIFR{i}{j}(k,1:end-1), GaussWin, 'same');
            end
            Temp_INIFR{i}{j}(:,1:end-1) = Temp_INIFR{i}{j}(:,1:end-1) - repmat(mean(Temp_INIFR{i}{j}(:,1:end-1), 2), 1, size(Temp_INIFR{i}{j}(:,1:end-1),2));
            
            clear NormTempINIFR;
            for k = 1:size(Temp_INIFR{i}{j},1),
                NormTempINIFR(k) = norm(Temp_INIFR{i}{j}(k,1:end-1));
            end
            
            TempCorr = [];
            for k = 1:size(LastINST,1),
                TempCorr = [TempCorr (LastINST(k,:)*Temp_INIFR{i}{j}(:,1:end-1)')./(NormLastINST(k)*NormTempINIFR)];
            end
            CorrToLastIN{i}(j) = mean(TempCorr(find(~isnan(TempCorr))));
            STDCorrToLastIN{i}(j) = std(TempCorr(find(~isnan(TempCorr))));
            
            TempCorr = [];
            for k = 1:size(FirstINST,1),
                TempCorr = [TempCorr (FirstINST(k,:)*Temp_INIFR{i}{j}(:,1:end-1)')./(NormFirstINST(k)*NormTempINIFR)];
            end
            CorrToFirstIN{i}(j) = mean(TempCorr(find(~isnan(TempCorr))));
            STDCorrToFirstIN{i}(j) = std(TempCorr(find(~isnan(TempCorr))));
        end
    end
end

MiddleINs = find(AllINPosition(:,2) == 1);
MiddleINST = AllINIFR(MiddleINs,1:end-1);
for i = 1:size(MiddleINST, 1),
    MiddleINST(i,:) = conv(MiddleINST(i,:), GaussWin, 'same');
end
MiddleINST = MiddleINST - repmat(mean(MiddleINST, 2), 1, size(MiddleINST,2));
NormMiddleINST = [];
for i = 1:size(MiddleINST,1),
    NormMiddleINST(i) = norm(MiddleINST(i,:));
end
TempCorr = [];
for i = 1:size(MiddleINST,1),
    TempCorr = [TempCorr (MiddleINST(i,:)*MiddleINST(i+1:end,:)')./(NormMiddleINST(i)*NormMiddleINST(i+1:end))];
end
INNeuralAnalysisResults.MiddleINs_Corr = [mean(TempCorr(find(~isnan(TempCorr)))) std(TempCorr(find(~isnan(TempCorr))))];


% Now to compare First INs with Last INs - pairwise correlation
TempCorr = [];
for i = 1:size(FirstINST,1),
    TempCorr = [TempCorr (FirstINST(i,:)*LastINST')./(NormFirstINST(i)*NormLastINST)];
end
INNeuralAnalysisResults.FirstINLastINCorr = [mean(TempCorr(find(~isnan(TempCorr)))) std(TempCorr(find(~isnan(TempCorr))))];

% Now to compare First INs with Middle INs - pairwise correlation
if (~isempty(MiddleINs))
    TempCorr = [];
    for i = 1:size(FirstINST,1),
        TempCorr = [TempCorr (FirstINST(i,:)*MiddleINST')./(NormFirstINST(i)*NormMiddleINST)];
    end
    INNeuralAnalysisResults.FirstINMiddleINCorr = [mean(TempCorr(find(~isnan(TempCorr)))) std(TempCorr(find(~isnan(TempCorr))))];
else
    INNeuralAnalysisResults.FirstINMiddleINCorr = NaN;
end

% Now to compare Last INs with Middle INs - pairwise correlation
if (~isempty(MiddleINs))
    TempCorr = [];
    for i = 1:size(MiddleINST,1),
        TempCorr = [TempCorr (MiddleINST(i,:)*LastINST')./(NormMiddleINST(i)*NormLastINST)];
    end
    INNeuralAnalysisResults.MiddleINLastINCorr = [mean(TempCorr(find(~isnan(TempCorr)))) std(TempCorr(find(~isnan(TempCorr))))];
else
    INNeuralAnalysisResults.MiddleINLastINCorr = NaN;
end

INNeuralAnalysisResults.INRaster = INRaster;
INNeuralAnalysisResults.INPST = INPST;
INNeuralAnalysisResults.INFR = INFR;
INNeuralAnalysisResults.INFR_IFR = INFR_IFR;

for i = 1:length(INNeuralAnalysisResults.INFR_IFR),
    if (~isempty(INNeuralAnalysisResults.INFR_IFR{i})),
        INNeuralAnalysisResults.NormINFR_IFR{i} = (INNeuralAnalysisResults.INFR_IFR{i} - mean(AllINFR_IFR))/std(AllINFR_IFR);
    end
end
INNeuralAnalysisResults.AllINPreMotorFR = AllINPreMotorFR;
INNeuralAnalysisResults.AllINPreMotorNumSpikes = AllINPreMotorNumSpikes;
INNeuralAnalysisResults.AllINIntPreMotorFR = AllINIntPreMotorFR;
INNeuralAnalysisResults.AllINIntLen = AllINIntLen;
INNeuralAnalysisResults.AllINIntPreMotorNumSpikes = AllINIntPreMotorNumSpikes;

INNeuralAnalysisResults.AllINSpikeTrain = AllINSpikeTrain;
INNeuralAnalysisResults.AllINPosition = AllINPosition;
INNeuralAnalysisResults.AllINIFR = AllINIFR;
INNeuralAnalysisResults.AllINFR = AllINFR;
INNeuralAnalysisResults.AllINFR_IFR = AllINFR_IFR;
INNeuralAnalysisResults.AllINPST = AllINPST;
INNeuralAnalysisResults.AllINFeats = AllINFeats;
% INNeuralAnalysisResults.DistanceToLast = pdist2(AllINFeats, mean(AllINFeats(find(AllINPosition(:,7) == -1),:)), 'mahalanobis', cov(AllINFeats(find(AllINPosition(:,7) == -1),:)));

if (exist('Corr', 'var'))
    INNeuralAnalysisResults.Corr = Corr;
    INNeuralAnalysisResults.StdCorr = StdCorr;
    INNeuralAnalysisResults.CorrToFirstIN = CorrToFirstIN;
    INNeuralAnalysisResults.STDCorrToFirstIN = STDCorrToFirstIN;
    INNeuralAnalysisResults.CorrToLastIN = CorrToLastIN;
    INNeuralAnalysisResults.STDCorrToLastIN = STDCorrToLastIN;
else
    INNeuralAnalysisResults.Corr = [];
    INNeuralAnalysisResults.StdCorr = [];
    INNeuralAnalysisResults.CorrToFirstIN = [];
    INNeuralAnalysisResults.STDCorrToFirstIN = [];
    INNeuralAnalysisResults.CorrToLastIN = [];
    INNeuralAnalysisResults.STDCorrToLastIN = [];
    INNeuralAnalysisResults.MeanLastINCorr = [];
    INNeuralAnalysisResults.MeanFirstINCorr = [];
end
INNeuralAnalysisResults.INFR = INFR;
INNeuralAnalysisResults.INFR_IFR = INFR_IFR;
INNeuralAnalysisResults.INRaster = INRaster;

INNeuralAnalysisResults.Edges = Edges;
INNeuralAnalysisResults.IFREdges = IFREdges;

INNeuralAnalysisResults.FirstINs = FirstINs;
INNeuralAnalysisResults.FirstINST = FirstINST;

INNeuralAnalysisResults.MiddleINs = MiddleINs;
INNeuralAnalysisResults.MiddleINST = MiddleINST;

INNeuralAnalysisResults.LastINs = LastINs;
INNeuralAnalysisResults.LastINST = LastINST;

% Now to compare correlation between mean patterns of activity for first,
% middle and last INs

INNeuralAnalysisResults.MeanFirstINLastINCorr = mean(FirstINST)*mean(LastINST)'./(norm(mean(FirstINST))*norm(mean(LastINST)));
INNeuralAnalysisResults.MeanFirstINMiddleINCorr = mean(FirstINST)*mean(MiddleINST)'./(norm(mean(FirstINST))*norm(mean(MiddleINST)));
INNeuralAnalysisResults.MeanLastINMiddleINCorr = mean(LastINST)*mean(MiddleINST)'./(norm(mean(LastINST))*norm(mean(MiddleINST)));

% Now to compare First INs with Middle INs - pairwise correlation

Temp_INIFR = INIFR;

MeanLastINST = mean(LastINST);
MeanFirstINST = mean(FirstINST);

for i = 1:MaxINs,
    for j = 1:length(Temp_INIFR{i}),
        if (size(Temp_INIFR{i}{j},1) > MinTrialNumber)
            FirstINCorr{i}(j) = 0;
            LastINCorr{i}(j) = 0;
            TempCorr = [];
            for k = 1:size(Temp_INIFR{i}{j},1),
                Temp_INIFR{i}{j}(k,1:end-1) = conv(Temp_INIFR{i}{j}(k,1:end-1), GaussWin, 'same');
            end
            Temp_INIFR{i}{j}(:,1:end-1) = Temp_INIFR{i}{j}(:,1:end-1) - repmat(mean(Temp_INIFR{i}{j}(:,1:end-1), 2), 1, size(Temp_INIFR{i}{j}(:,1:end-1),2));
            for k = 1:size(Temp_INIFR{i}{j},1),
                Norm_Temp_INIFR{i}{j}(k,1) = norm(Temp_INIFR{i}{j}(k,1:end-1));
            end
            
            TempCorr = (Temp_INIFR{i}{j}(:,1:end-1)*MeanFirstINST')./(Norm_Temp_INIFR{i}{j}*norm(MeanFirstINST));
            FirstINCorr{i}(j) = mean(TempCorr);
            FirstINStdCorr{i}(j) = std(TempCorr);
            
            TempCorr = (Temp_INIFR{i}{j}(:,1:end-1)*MeanLastINST')./(Norm_Temp_INIFR{i}{j}*norm(MeanLastINST));
            LastINCorr{i}(j) = mean(TempCorr);
            LastINStdCorr{i}(j) = std(TempCorr);
            
            INNeuralAnalysisResults.MeanLastINCorr{i}(j) = mean(LastINST)*mean(Temp_INIFR{i}{j}(:,1:end-1))'./(norm(mean(LastINST))*norm(mean(Temp_INIFR{i}{j}(:,1:end-1))));
            INNeuralAnalysisResults.MeanFirstINCorr{i}(j) = mean(FirstINST)*mean(Temp_INIFR{i}{j}(:,1:end-1))'./(norm(mean(FirstINST))*norm(mean(Temp_INIFR{i}{j}(:,1:end-1))));
            
        end
    end
end
if (~exist('FirstINCorr', 'var'))
    FirstINCorr = [];
    LastINCorr = [];
    FirstINStdCorr = [];
    LastINStdCorr = [];
end
    
INNeuralAnalysisResults.FirstINCorr = FirstINCorr;
INNeuralAnalysisResults.FirstINStdCorr = FirstINStdCorr;
INNeuralAnalysisResults.LastINCorr = LastINCorr;
INNeuralAnalysisResults.LastINStdCorr = LastINStdCorr;

if (strfind(PlotOption, 'on'))

    if (~isempty(Indices))
        PlotHGap = 0.03;
        PlotWGap = 0.075;

        PlotHt = 0.92/length(Indices) - PlotHGap;

        PlotWidth = 0.96/max(Indices) - PlotWGap;

        RasterFig = figure;
        set(gcf, 'Color', 'w');

        FRFig = figure;
        set(gcf, 'Color', 'w');

        CorrFig = figure;
        set(gcf, 'Color', 'w');

        RowNo = 0;
    
        for i = Indices,
            for j = 1:length(INRaster{i})-1,
                figure(RasterFig);
                subplot('Position', [(1 - (length(INRaster{i}) - j)*(PlotWidth + PlotWGap)) (1 - (RowNo+1)*(PlotHt + PlotHGap)) PlotWidth PlotHt]);
                %subplot(length(Indices), max(Indices), RowNo + j + max(Indices) - length(INRaster{i}) + 1);
                if (strfind(PlotMinINs, 'no'))
                    MinTrials = size(INPST{i}{j},1);
                end
                hold on;
                if (strfind(PlotPSTH, 'on'))
                    plot(Edges(1:end-1), mean(INPST{i}{j}(:,1:end-1))/max([max(AllINPST(:,1:end-1))]) * (MinTrials), 'Color', 'r', 'LineWidth', 1.5);
                end
                %MeanPST = mean(INPST{i}{j}(:,1:end-1))/max([max(AllINPST(:,1:end-1))]) * MinTrials;
                %STDPST = std(INPST{i}{j}(:,1:end-1))/max([max(AllINPST(:,1:end-1))]) * MinTrials;
                %SEMPST = (std(INPST{i}{j}(:,1:end-1))/sqrt(size(INPST{i}{j}(:,1:end-1),1)))/max([max(AllINPST(:,1:end-1))]) * MinTrials;
                %fill([Edges(1:end-1)+0.00001 fliplr(Edges(1:end-1)-0.00001)], [MeanPST+SEMPST fliplr(MeanPST-SEMPST)], 'r', 'EdgeColor', 'none', 'FaceColor', 'r', 'FaceAlpha', 0.75);
                if (~isempty(INRaster{i}{j}))
                    PlotRaster(INRaster{i}{j}, 'k', 0.25, 0, MinTrials);
                end
                fill([zeros(1,MinTrials) flipud(SyllOffsetRaster{i}{j}(1:MinTrials,1))'], [linspace(0.5, MinTrials+0.5, MinTrials) linspace(MinTrials+0.5, 0.5, MinTrials)], 'k', 'EdgeColor', 'none', 'FaceColor', [0.65 0.65 0.65], 'FaceAlpha', 0.5);
                %plot([0 0], [-0.5 MinTrials+0.5], 'k--');
                %plot(SyllOffsetRaster{i}{j}(1:MinTrials,1), ones(MinTrials, 1), 'ro', 'MarkerSize', 3);
                %fill([(zeros(1,MinTrials)-0.045) flipud(SyllOffsetRaster{i}{j}(1:MinTrials,1) - 0.045)'], [linspace(0.5, MinTrials+0.5, MinTrials) linspace(MinTrials+0.5, 0.5, MinTrials)], 'k', 'EdgeColor', 'none', 'FaceColor', [1 0.75 0.75], 'FaceAlpha', 0.5);
                %fill([(SyllOffsetRaster{i}{j}(1:MinTrials,1) - 0.045)' flipud(NextSyllOnsetRaster{i}{j}(1:MinTrials,1) - 0.045)'], [linspace(0.5, MinTrials+0.5, MinTrials) linspace(MinTrials+0.5, 0.5, MinTrials)], 'k', 'EdgeColor', 'none', 'FaceColor', [0.75 0.75 1], 'FaceAlpha', 0.5);
                %plot(PrevSyllOffsetRaster{i}{j}(1:MinTrials,1), PrevSyllOffsetRaster{i}{j}(1:MinTrials,2), 'ro', 'MarkerSize', 4, 'MarkerFace', 'r');
                %plot(SyllOffsetRaster{i}{j}(1:MinTrials,1), SyllOffsetRaster{i}{j}(1:MinTrials,2), 'ro', 'MarkerSize', 2);
                %plot(NextSyllOnsetRaster{i}{j}(1:MinTrials,1), NextSyllOnsetRaster{i}{j}(1:MinTrials,2), 'ro', 'MarkerSize', 4, 'MarkerFace', 'r');
                axis([Edges(1) Edges(end) 0.5 MinTrials+0.5]);
                %plot([0 0], [0.5 MinTrials+0.5], 'g--', 'LineWidth', 2);
                if (i ~= Indices(end))
                    set(gca, 'XTick', [-0.05 0 0.05], 'XTickLabel', []);
                else
                    set(gca, 'XTick', [-0.05 0 0.05], 'XTickLabel', [-50 0 50]);
                end
                
                set(gca, 'Box', 'off');
                set(gca, 'TickLen', [0.05 0.025]);
                set(gca, 'FontSize', 18, 'FontName', 'Times New Roman');
                % set(gca, 'TickDir', 'out');
                set(gca, 'YTick', []);
                set(gca, 'YAxisLocation', 'right');
                set(gca, 'YColor', 'w');

%                 if (i == Indices(end))
%                     % X-axis with tick marks
%                     plot([-0.05 0.05], [-0.5 -0.5], 'k');
%                     for TickMark = -0.05:0.05:0.05,
%                         plot([TickMark TickMark], [-0.5 0], 'k');
%                     end
%                 end
                if (strfind(PlotPSTH, 'on'))
                    if (j == 1)
                        plot([-0.05 -0.05], [(MinTrials+0.5) ((MinTrials+0.5) - (100/max([max(AllINPST(:,1:end-1))]) * MinTrials))], 'r');
                    end
                end
            end
            
            % Now to plot the firing rates
            figure(FRFig);
            hold on;
            set(gcf, 'Position', [920 277 150 420]);
            FRAxis(i-min(Indices)+1) = subplot('Position', [0.22 (1 - (RowNo+1)*(PlotHt + PlotHGap)) 0.76 PlotHt]);
            %errorbar([-i:1:-1], mean(INFR_IFR{i}(:,1:end-1)), std(INFR_IFR{i}(:,1:end-1)), 'ks-'); % firing rates using mean IFR
            errorbar([-i:1:-1], mean(INFR{i}(:,1:end-1)), std(INFR{i}(:,1:end-1)), 'ks-'); % firing rates using mean spike count
            axis tight;
            temp = axis;
            axis([(-max(Indices)-0.25) -0.75 0 temp(4)*1.1]);
            FinalAxis(i-min(Indices)+1,:) = axis;
            set(gca, 'Box', 'off');
            set(gca, 'XTick', []);
            set(gca, 'FontSize', 8, 'FontName', 'Times New Roman');
            %set(gca, 'XTick', [-i:1:-1]);
            set(gca, 'TickLength', [0.015 0.025]);
            
            % Now to plot the correlations
            figure(CorrFig);
            set(gcf, 'Position', [1070 277 150 420]);
            subplot('Position', [0.22 (1 - (RowNo+1)*(PlotHt + PlotHGap)) 0.76 PlotHt]);
            %errorbar([-i:1:-1], CorrToFirstIN{i}(1:end-1), STDCorrToFirstIN{i}(1:end-1), 'rs-');
            %plot([-i:1:-1], INNeuralAnalysisResults.MeanFirstINCorr{i}(1:end-1), 'rs-');
            hold on;
            %errorbar([-i:1:-1], CorrToLastIN{i}(1:end-1), STDCorrToLastIN{i}(1:end-1), 'bs-');
            plot([-i:1:-1], INNeuralAnalysisResults.MeanLastINCorr{i}(1:end-1), 'ks-');
            % errorbar([-i:1:-1], LastINCorr{i}(1:end-1), LastINStdCorr{i}(1:end-1), 'rd-');

            axis([(-max(Indices)-0.25) -0.75 -1 1]);
            plot([(-max(Indices)-0.25) -0.75], [0 0], 'k--');
            
            set(gca, 'Box', 'off');
            set(gca, 'XTick', []);
            set(gca, 'FontSize', 8, 'FontName', 'Times New Roman');
            set(gca, 'YTick', [-1 0 1]);
            %set(gca, 'XTick', [-i:1:-1]);
            set(gca, 'TickLength', [0.015 0.025]);
            RowNo = RowNo + 1;
        end
    end

    FinalFRAxis = [min(FinalAxis(:,1)) max(FinalAxis(:,2)) 0 max(FinalAxis(:,4))];
    for i = 1:length(FRAxis),
        if (FRAxis(i) > 0)
            axes(FRAxis(i));
            axis(FinalFRAxis);
            set(gca, 'YTick', [0:50:FinalFRAxis(end)], 'YTickLabel', []);
            set(gca, 'TickLen', [0.025 0.025]);
        end
    end
    
    disp(['Firing rate scale bar is 100 Hz']);
end

INNeuralAnalysisResults.LastvsOthersIndex = (mean(AllINFR(find(AllINPosition(:,end) == -1))) - mean(AllINFR(find(AllINPosition(:,end) ~= -1))))/(mean(AllINFR(find(AllINPosition(:,end) == -1))) + mean(AllINFR(find(AllINPosition(:,end) ~= -1))));

INNeuralAnalysisResults.LastvsOthersFR(1) = mean(AllINFR(find(AllINPosition(:,end) == -1)));
INNeuralAnalysisResults.LastvsOthersFR(2) = std(AllINFR(find(AllINPosition(:,end) == -1)))/sqrt(length(find(AllINPosition(:,end) == -1)));
INNeuralAnalysisResults.LastvsOthersFR(3) = mean(AllINFR(find(AllINPosition(:,end) ~= -1)));
INNeuralAnalysisResults.LastvsOthersFR(4) = std(AllINFR(find(AllINPosition(:,end) ~= -1)))/sqrt(length(find(AllINPosition(:,end) == -1)));

disp('Finished Analysis');