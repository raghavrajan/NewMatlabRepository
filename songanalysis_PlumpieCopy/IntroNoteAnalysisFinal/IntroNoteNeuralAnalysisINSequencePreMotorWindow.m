function [INNeuralAnalysisResults] = IntroNoteNeuralAnalysisINSequencePreMotorWindow(Neural_INR, PlotOption, MinTrialNumber, PlotPSTH)

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

PreMotorOffset = 0.045;

Width = 0.004;
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
BinSize = 0.001;
LastGapDur = [];
for i = 1:length(Neural_INR.BoutDetails),
    Matches = strfind(Neural_INR.BoutDetails(i).labels, 'ia');
    for j = 1:length(Matches),
        LastGapDur = [LastGapDur; (Neural_INR.BoutDetails(i).onsets(Matches(j)+1) - Neural_INR.BoutDetails(i).offsets(Matches(j)))];
    end
end
MedianLastGapDur = median(LastGapDur);
LastGapDur = min(LastGapDur);
disp(['The minimum last gap duration is ', num2str(LastGapDur)]);

INDur = [];
INGap = [];
Index = 0;
for i = 1:length(Neural_INR.BoutDetails),
    Index = Index + 1;
    if (Neural_INR.NoofINs(i) > 0)
        Neural_INR.INDurs{Index} = Neural_INR.BoutDetails(i).offsets(Neural_INR.INs{i}) - Neural_INR.BoutDetails(i).onsets(Neural_INR.INs{i});
        Neural_INR.INGaps{Index} = Neural_INR.BoutDetails(i).onsets(Neural_INR.INs{i}+1) - Neural_INR.BoutDetails(i).offsets(Neural_INR.INs{i});
        INDur = [INDur; Neural_INR.INDurs{i}];
        INGap = [INGap; Neural_INR.INGaps{i}];
    else
        Neural_INR.INDurs{Index} = [];
        Neural_INR.INGaps{Index} = [];
    end
end

for i = 1:size(Neural_INR.WithinBoutNoofINs,1),
    Index = Index + 1;
    if (Neural_INR.WithinBoutNoofINs(i,1) > 0)
        Neural_INR.WithinBoutINDurs{Index} = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).offsets(Neural_INR.WithinBoutINs{i}) - Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(Neural_INR.WithinBoutINs{i});
        Neural_INR.WithinBoutINGaps{Index} = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(Neural_INR.WithinBoutINs{i}+1) - Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).offsets(Neural_INR.WithinBoutINs{i});
        INDur = [INDur; Neural_INR.WithinBoutINDurs{i}];
        INGap = [INGap; Neural_INR.WithinBoutINGaps{i}];
    else
        Neural_INR.WithinBoutINDurs{Index} = [];
        Neural_INR.WithinBoutINDurs{Index} = [];
    end
end

INPSTEdges = -PreMotorOffset:BinSize:(median(INDur) - PreMotorOffset);
GapPSTEdges = (median(INDur) - PreMotorOffset):BinSize:(median(INDur) - PreMotorOffset + LastGapDur);

disp(['The median IN duration is ', num2str(median(INDur)), ' and the median gap duration is ', num2str(median(INGap))]);

MaxINs = max([max(Neural_INR.NoofINs) max(Neural_INR.WithinBoutNoofINs(:,1))]);

for i = 1:MaxINs,
    for j = 1:i,
        INRaster{i}{j} = [];
        INPST{i}{j} = [];
        GapPST{i}{j} = [];
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
            
            for k = 1:i,
                AllINIndex = AllINIndex + 1;
                Edges = ((StartTimes(k) - PreMotorOffset):BinSize:(StartTimes(k+1) - PreMotorOffset)) - StartTimes(k);
                
                INIFRStartIndex = find(Neural_INR.BoutDetails(j).IFR(1,:) <= (StartTimes(k) - PreMotorOffset), 1, 'last');
                INIFREndIndex = find(Neural_INR.BoutDetails(j).IFR(1,:) <= (EndTimes(k) - PreMotorOffset), 1, 'last');
                INIFR{i}{k}{TrialNo} = interp1(Neural_INR.BoutDetails(j).IFR(1, INIFRStartIndex:INIFREndIndex), Neural_INR.BoutDetails(j).IFR(2, INIFRStartIndex:INIFREndIndex), ((StartTimes(k) - PreMotorOffset):1/IFRFs:(EndTimes(k) - PreMotorOffset)));

                GapIFRStartIndex = find(Neural_INR.BoutDetails(j).IFR(1,:) <= (EndTimes(k) - PreMotorOffset), 1, 'last');
                GapIFREndIndex = find(Neural_INR.BoutDetails(j).IFR(1,:) <= (StartTimes(k+1) - PreMotorOffset), 1, 'last');
                GapIFR{i}{k}{TrialNo} = interp1(Neural_INR.BoutDetails(j).IFR(1, INIFRStartIndex:INIFREndIndex), Neural_INR.BoutDetails(j).IFR(2, INIFRStartIndex:INIFREndIndex), ((StartTimes(k) - PreMotorOffset):1/IFRFs:(EndTimes(k) - PreMotorOffset)));
                
                INPreMotorSpikeIndices = find((Neural_INR.BoutDetails(j).SpikeTimes >= (StartTimes(k) - PreMotorOffset)) & (Neural_INR.BoutDetails(j).SpikeTimes < (EndTimes(k) - PreMotorOffset)));
                
                GapPreMotorSpikeIndices = find((Neural_INR.BoutDetails(j).SpikeTimes >= (EndTimes(k) - PreMotorOffset)) & (Neural_INR.BoutDetails(j).SpikeTimes < (EndTimes(k) - PreMotorOffset + LastGapDur)));
                FullGapPreMotorSpikeIndices = find((Neural_INR.BoutDetails(j).SpikeTimes >= (EndTimes(k) - PreMotorOffset)) & (Neural_INR.BoutDetails(j).SpikeTimes < (StartTimes(k+1) - PreMotorOffset)));
                
                if (~isempty(INPreMotorSpikeIndices))
                    SpikeTimes = ((Neural_INR.BoutDetails(j).SpikeTimes(INPreMotorSpikeIndices) - (StartTimes(k) - PreMotorOffset)) * median(INDur)/(EndTimes(k) - StartTimes(k))) - PreMotorOffset;
                    INRaster{i}{k} = [INRaster{i}{k}; [SpikeTimes ones(length(INPreMotorSpikeIndices),1)*TrialNo]];
                    AllINPreMotorSpikeTrain{AllINIndex,:} = SpikeTimes;
                    INPST{i}{k}(TrialNo,:) = histc(SpikeTimes, INPSTEdges);
                else
                    INPST{i}{k}(TrialNo,:) = zeros(1, length(INPSTEdges));
                    AllINPreMotorSpikeTrain{AllINIndex,:} = [];
                end
                
                if (~isempty(GapPreMotorSpikeIndices))
                    SpikeTimes = (Neural_INR.BoutDetails(j).SpikeTimes(GapPreMotorSpikeIndices) - EndTimes(k) + median(INDur));
                    INRaster{i}{k} = [INRaster{i}{k}; [SpikeTimes ones(length(SpikeTimes),1)*TrialNo]];
                    AllGapPreMotorSpikeTrain{AllINIndex,:} = SpikeTimes;
                    GapPST{i}{k}(TrialNo,:) = histc(SpikeTimes, GapPSTEdges);
                else
                    GapPST{i}{k}(TrialNo,:) = zeros(1, length(GapPSTEdges));
                    AllGapPreMotorSpikeTrain{AllINIndex,:} = [];
                end
                
                if (~isempty(FullGapPreMotorSpikeIndices))
                    % FullGapWarpedSpikeTimes = ((Neural_INR.BoutDetails(j).SpikeTimes(FullGapPreMotorSpikeIndices) - (EndTimes(k) - PreMotorOffset)) * median(INGap)/(StartTimes(k+1) - EndTimes(k))) + median(INDur) - PreMotorOffset;
                    FullGapWarpedSpikeTimes = Neural_INR.BoutDetails(j).SpikeTimes(FullGapPreMotorSpikeIndices) - EndTimes(k) + median(INDur);
                    % INRaster{i}{k} = [INRaster{i}{k}; [FullGapWarpedSpikeTimes ones(length(FullGapWarpedSpikeTimes),1)*TrialNo]];
                    AllFullGapPreMotorSpikeTrain{AllINIndex,:} = FullGapWarpedSpikeTimes;
                else
                    AllFullGapPreMotorSpikeTrain{AllINIndex,:} = [];
                end
               
                INPreMotorFR{i}(TrialNo, k) = length(INPreMotorSpikeIndices)/(EndTimes(k) - StartTimes(k));
                INGapPreMotorFR{i}(TrialNo, k) = length(GapPreMotorSpikeIndices)/LastGapDur;
                INFullGapPreMotorFR{i}(TrialNo, k) = length(FullGapPreMotorSpikeIndices)/(StartTimes(k+1) - EndTimes(k));
                PrevSyllOffsetRaster{i}{k} = [PrevSyllOffsetRaster{i}{k}; [(PrevSyllOffsetTimes(k)-StartTimes(k)) TrialNo]];
                SyllOffsetRaster{i}{k} = [SyllOffsetRaster{i}{k}; [(median(INDur)) TrialNo]];
                NextSyllOnsetRaster{i}{k} = [NextSyllOnsetRaster{i}{k}; [(StartTimes(k+1) - EndTimes(k) - PreMotorOffset + median(INDur)) TrialNo]];

                AllINPreMotorFR(AllINIndex,:) = length(INPreMotorSpikeIndices)/(EndTimes(k) - StartTimes(k));
                AllINPreMotorNumSpikes(AllINIndex,:) = length(INPreMotorSpikeIndices);
                
                AllINGapPreMotorFR(AllINIndex,:) = length(GapPreMotorSpikeIndices)/LastGapDur;
                AllINGapPreMotorNumSpikes(AllINIndex,:) = length(GapPreMotorSpikeIndices);
                
                AllINFullGapLen(AllINIndex,:) = StartTimes(k+1) - EndTimes(k);
                AllINFullGapPreMotorFR(AllINIndex,:) = length(FullGapPreMotorSpikeIndices)/(StartTimes(k+1) - EndTimes(k));
                AllINFullGapPreMotorNumSpikes(AllINIndex,:) = length(FullGapPreMotorSpikeIndices);
               
                AllINPST(AllINIndex,:) = [INPST{i}{k}(TrialNo,:) GapPST{i}{k}(TrialNo,:)];
                AllINFeats(AllINIndex,:) = Neural_INR.BoutDetails(j).Feats(Neural_INR.INs{j}(k),1:4);
                
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
            
            for k = 1:i,
                AllINIndex = AllINIndex + 1;
                Edges = ((StartTimes(k) - PreMotorOffset):BinSize:(StartTimes(k+1) - PreMotorOffset)) - StartTimes(k);
                
                INIFRStartIndex = find(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).IFR(1,:) <= (StartTimes(k) - PreMotorOffset), 1, 'last');
                INIFREndIndex = find(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).IFR(1,:) <= (EndTimes(k) - PreMotorOffset), 1, 'last');
                INIFR{i}{k}{TrialNo} = interp1(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).IFR(1, INIFRStartIndex:INIFREndIndex), Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).IFR(2, INIFRStartIndex:INIFREndIndex), ((StartTimes(k) - PreMotorOffset):1/IFRFs:(EndTimes(k) - PreMotorOffset)));

                GapIFRStartIndex = find(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).IFR(1,:) <= (EndTimes(k) - PreMotorOffset), 1, 'last');
                GapIFREndIndex = find(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).IFR(1,:) <= (StartTimes(k+1) - PreMotorOffset), 1, 'last');
                GapIFR{i}{k}{TrialNo} = interp1(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).IFR(1, INIFRStartIndex:INIFREndIndex), Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).IFR(2, INIFRStartIndex:INIFREndIndex), ((StartTimes(k) - PreMotorOffset):1/IFRFs:(EndTimes(k) - PreMotorOffset)));
                
                INPreMotorSpikeIndices = find((Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).SpikeTimes >= (StartTimes(k) - PreMotorOffset)) & (Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).SpikeTimes < (EndTimes(k) - PreMotorOffset)));
                
                GapPreMotorSpikeIndices = find((Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).SpikeTimes >= (EndTimes(k) - PreMotorOffset)) & (Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).SpikeTimes < (EndTimes(k) - PreMotorOffset + LastGapDur)));
                FullGapPreMotorSpikeIndices = find((Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).SpikeTimes >= (EndTimes(k) - PreMotorOffset)) & (Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).SpikeTimes < (StartTimes(k+1) - PreMotorOffset)));
                
                if (~isempty(INPreMotorSpikeIndices))
                    SpikeTimes = ((Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).SpikeTimes(INPreMotorSpikeIndices) - (StartTimes(k) - PreMotorOffset)) * median(INDur)/(EndTimes(k) - StartTimes(k))) - PreMotorOffset;
                    INRaster{i}{k} = [INRaster{i}{k}; [SpikeTimes ones(length(INPreMotorSpikeIndices),1)*TrialNo]];
                    AllINPreMotorSpikeTrain{AllINIndex,:} = SpikeTimes;
                    INPST{i}{k}(TrialNo,:) = histc(SpikeTimes, INPSTEdges);
                else
                    INPST{i}{k}(TrialNo,:) = zeros(1, length(INPSTEdges));
                    AllINPreMotorSpikeTrain{AllINIndex,:} = [];
                end
                
                if (~isempty(GapPreMotorSpikeIndices))
                    SpikeTimes = (Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).SpikeTimes(GapPreMotorSpikeIndices) - EndTimes(k) + median(INDur));
                    INRaster{i}{k} = [INRaster{i}{k}; [SpikeTimes ones(length(SpikeTimes),1)*TrialNo]];
                    AllGapPreMotorSpikeTrain{AllINIndex,:} = SpikeTimes;
                    GapPST{i}{k}(TrialNo,:) = histc(SpikeTimes, GapPSTEdges);
                else
                    GapPST{i}{k}(TrialNo,:) = zeros(1, length(GapPSTEdges));
                    AllGapPreMotorSpikeTrain{AllINIndex,:} = [];
                end
                
                if (~isempty(FullGapPreMotorSpikeIndices))
                    % FullGapWarpedSpikeTimes = ((Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).SpikeTimes(FullGapPreMotorSpikeIndices) - (EndTimes(k) - PreMotorOffset)) * median(INGap)/(StartTimes(k+1) - EndTimes(k))) + median(INDur) - PreMotorOffset;
                    FullGapWarpedSpikeTimes = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).SpikeTimes(FullGapPreMotorSpikeIndices) - EndTimes(k) + median(INDur);
                    % INRaster{i}{k} = [INRaster{i}{k}; [FullGapWarpedSpikeTimes ones(length(FullGapWarpedSpikeTimes),1)*TrialNo]];
                    AllFullGapPreMotorSpikeTrain{AllINIndex,:} = FullGapWarpedSpikeTimes;
                else
                    AllFullGapPreMotorSpikeTrain{AllINIndex,:} = [];
                end
                
                INPreMotorFR{i}(TrialNo, k) = length(INPreMotorSpikeIndices)/(EndTimes(k) - StartTimes(k));
                INGapPreMotorFR{i}(TrialNo, k) = length(GapPreMotorSpikeIndices)/LastGapDur;
                INFullGapPreMotorFR{i}(TrialNo, k) = length(FullGapPreMotorSpikeIndices)/(StartTimes(k+1) - EndTimes(k));
                PrevSyllOffsetRaster{i}{k} = [PrevSyllOffsetRaster{i}{k}; [(PrevSyllOffsetTimes(k)-StartTimes(k)) TrialNo]];
                SyllOffsetRaster{i}{k} = [SyllOffsetRaster{i}{k}; [(median(INDur)) TrialNo]];
                NextSyllOnsetRaster{i}{k} = [NextSyllOnsetRaster{i}{k}; [(StartTimes(k+1) - EndTimes(k) - PreMotorOffset + median(INDur)) TrialNo]];

                AllINPreMotorFR(AllINIndex,:) = length(INPreMotorSpikeIndices)/(EndTimes(k) - StartTimes(k));
                AllINPreMotorNumSpikes(AllINIndex,:) = length(INPreMotorSpikeIndices);
                
                AllINGapPreMotorFR(AllINIndex,:) = length(GapPreMotorSpikeIndices)/LastGapDur;
                AllINGapPreMotorNumSpikes(AllINIndex,:) = length(GapPreMotorSpikeIndices);
                
                AllINFullGapLen(AllINIndex,:) = StartTimes(k+1) - EndTimes(k);
                AllINFullGapPreMotorFR(AllINIndex,:) = length(FullGapPreMotorSpikeIndices)/(StartTimes(k+1) - EndTimes(k));
                AllINFullGapPreMotorNumSpikes(AllINIndex,:) = length(FullGapPreMotorSpikeIndices);

                AllINPST(AllINIndex,:) = [INPST{i}{k}(TrialNo,:) GapPST{i}{k}(TrialNo,:)];
                AllINFeats(AllINIndex,:) = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).Feats(Neural_INR.WithinBoutINs{j}(k),1:4);
                
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

INNeuralAnalysisResults.AllINPreMotorFR = AllINPreMotorFR;
INNeuralAnalysisResults.AllINPreMotorNumSpikes = AllINPreMotorNumSpikes;
INNeuralAnalysisResults.AllINGapPreMotorFR = AllINGapPreMotorFR;
INNeuralAnalysisResults.AllINGapPreMotorNumSpikes = AllINGapPreMotorNumSpikes;
INNeuralAnalysisResults.AllINFullGapLen = AllINFullGapLen;
INNeuralAnalysisResults.AllINFullGapPreMotorFR = AllINFullGapPreMotorFR;
INNeuralAnalysisResults.AllINFullGapPreMotorNumSpikes = AllINFullGapPreMotorNumSpikes;
INNeuralAnalysisResults.AllINPreMotorSpikeTrain = AllINPreMotorSpikeTrain;
INNeuralAnalysisResults.AllGapPreMotorSpikeTrain = AllGapPreMotorSpikeTrain;
INNeuralAnalysisResults.AllFullGapPreMotorSpikeTrain = AllFullGapPreMotorSpikeTrain;

INNeuralAnalysisResults.INPreMotorFR = INPreMotorFR;
INNeuralAnalysisResults.INGapPreMotorFR = INGapPreMotorFR;
INNeuralAnalysisResults.INFullGapPreMotorFR = INFullGapPreMotorFR;
                
INNeuralAnalysisResults.AllINPosition = AllINPosition;
INNeuralAnalysisResults.AllINFeats = AllINFeats;

INNeuralAnalysisResults.INRaster = INRaster;

for i = min(AllINPosition(:,end)):1:max(AllINPosition(:,end)),
    TempIndices = find(AllINPosition(:,end) == i);
    if (length(TempIndices) >= MinTrials)
        break;
    end
end
MinPosition = i;

INNeuralAnalysisResults.ValidIndices = [];
for i = MinPosition:1:max(AllINPosition(:,end)),
    INNeuralAnalysisResults.ValidIndices = [INNeuralAnalysisResults.ValidIndices; find(AllINPosition(:,end) == i)];
end
    
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

   %     CorrFig = figure;
   %     set(gcf, 'Color', 'w');

        RowNo = 0;
    
        for i = Indices,
            for j = 1:length(INRaster{i}),
                figure(RasterFig);
                subplot('Position', [(1 - (length(INRaster{i}) + 1 - j)*(PlotWidth + PlotWGap)) (1 - (RowNo+1)*(PlotHt + PlotHGap)) PlotWidth PlotHt]);
                hold on;
                if (~isempty(INRaster{i}{j}))
                    PlotRaster(INRaster{i}{j}, 'k', 0.25, 0, MinTrials);
                end
                % plot([0 0], [-0.5 MinTrials+0.5], 'k');
%                plot(SyllOffsetRaster{i}{j}(1:MinTrials,1), ones(MinTrials, 1), 'ro', 'MarkerSize', 3);
                fill([(zeros(1,MinTrials) - PreMotorOffset) flipud(ones(1,MinTrials)*median(INDur) - PreMotorOffset)], [linspace(0.5, MinTrials+0.5, MinTrials) linspace(MinTrials+0.5, 0.5, MinTrials)], 'k', 'EdgeColor', 'none', 'FaceColor', [1 0.75 0.75], 'FaceAlpha', 0.5);
                fill([(ones(1,MinTrials)*median(INDur) - PreMotorOffset) flipud(ones(1, MinTrials)*median(INDur) - PreMotorOffset + LastGapDur)], [linspace(0.5, MinTrials+0.5, MinTrials) linspace(MinTrials+0.5, 0.5, MinTrials)], 'k', 'EdgeColor', 'none', 'FaceColor', [0.75 0.75 1], 'FaceAlpha', 0.5);                % fill([(ones(1,MinTrials)*median(INDur) - PreMotorOffset) flipud(NextSyllOnsetRaster{i}{j}(1:MinTrials,1))'], [linspace(0.5, MinTrials+0.5, MinTrials) linspace(MinTrials+0.5, 0.5, MinTrials)], 'k', 'EdgeColor', 'none', 'FaceColor', [0.75 0.75 1], 'FaceAlpha', 0.5);
                % fill([(ones(1,MinTrials)*median(INDur) - PreMotorOffset) flipud(ones(1,MinTrials)*median(INGap) - PreMotorOffset + median(INDur))], [linspace(0.5, MinTrials+0.5, MinTrials) linspace(MinTrials+0.5, 0.5, MinTrials)], 'k', 'EdgeColor', 'none', 'FaceColor', [0.75 0.75 1], 'FaceAlpha', 0.5);
                if (strfind(PlotPSTH, 'on'))
                    plot(INPSTEdges(:,1:end-1) + BinSize/2, mean(INPST{i}{j}(:,1:end-1))/max([max(AllINPST(:,1:end-1))]) * (MinTrials), 'Color', 'r', 'LineWidth', 1.5);
                    plot(GapPSTEdges(:,1:end-1) + BinSize/2, mean(GapPST{i}{j}(:,1:end-1))/max([max(AllINPST(:,1:end-1))]) * (MinTrials), 'Color', 'b', 'LineWidth', 1.5);
                end
                %plot(PrevSyllOffsetRaster{i}{j}(1:MinTrials,1), PrevSyllOffsetRaster{i}{j}(1:MinTrials,2), 'ro', 'MarkerSize', 4, 'MarkerFace', 'r');
                %plot(SyllOffsetRaster{i}{j}(1:MinTrials,1), SyllOffsetRaster{i}{j}(1:MinTrials,2), 'ro', 'MarkerSize', 2);
                %plot(NextSyllOnsetRaster{i}{j}(1:MinTrials,1), NextSyllOnsetRaster{i}{j}(1:MinTrials,2), 'ro', 'MarkerSize', 4, 'MarkerFace', 'r');
                axis([-PreMotorOffset (median(INDur) - PreMotorOffset + LastGapDur) -0.5 MinTrials+0.5]);
                temp = axis;
                if (i ~= Indices(end))
                    set(gca, 'XTick', [temp(1) 0], 'XTickLabel', []);
                else
                    set(gca, 'XTick', [temp(1) 0], 'XTickLabel', [temp(1)*1000 0]);
                end
                set(gca, 'Box', 'off');
                set(gca, 'TickLen', [0.05 0.025]);
                set(gca, 'FontSize', 8, 'FontName', 'Times New Roman');
                % set(gca, 'TickDir', 'out');
                set(gca, 'YTick', []);
                set(gca, 'YAxisLocation', 'right');
                set(gca, 'YColor', 'w');

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
            hold on;
            errorbar([-i:1:-1], mean(INPreMotorFR{i}(:,1:end)), std(INPreMotorFR{i}(:,1:end)), 'rs-'); % firing rates using mean spike count in the window premotor to the IN
            errorbar([-i:1:-1], mean(INGapPreMotorFR{i}(:,1:end)), std(INGapPreMotorFR{i}(:,1:end)), 'bs-'); % firing rates using mean spike count in the window premotor to the interval between two INs
            
            axis tight;
            temp = axis;
            axis([(-max(Indices)-0.25) -0.75 0 temp(4)*1.1]);
            FinalAxis(i-min(Indices)+1,:) = axis;
            set(gca, 'Box', 'off');
            set(gca, 'XTick', []);
            set(gca, 'FontSize', 8, 'FontName', 'Times New Roman');
            %set(gca, 'YTickLabel', []);
            %set(gca, 'XTick', [-i:1:-1]);
            set(gca, 'TickLength', [0.015 0.025]);
            
            % Now to plot the correlations
%             figure(CorrFig);
%             set(gcf, 'Position', [1070 277 150 420]);
%             subplot('Position', [0.22 (1 - (RowNo+1)*(PlotHt + PlotHGap)) 0.76 PlotHt]);
%             %errorbar([-i:1:-1], CorrToFirstIN{i}(1:end-1), STDCorrToFirstIN{i}(1:end-1), 'rs-');
%             %plot([-i:1:-1], INNeuralAnalysisResults.MeanFirstINCorr{i}(1:end-1), 'rs-');
%             hold on;
%             %errorbar([-i:1:-1], CorrToLastIN{i}(1:end-1), STDCorrToLastIN{i}(1:end-1), 'bs-');
%             plot([-i:1:-1], INNeuralAnalysisResults.MeanLastINCorr{i}(1:end-1), 'ks-');
%             % errorbar([-i:1:-1], LastINCorr{i}(1:end-1), LastINStdCorr{i}(1:end-1), 'rd-');
% 
%             axis([(-max(Indices)-0.25) -0.75 -1 1]);
%             plot([(-max(Indices)-0.25) -0.75], [0 0], 'k--');
%             
%             set(gca, 'Box', 'off');
%             set(gca, 'XTick', []);
%             set(gca, 'FontSize', 8, 'FontName', 'Times New Roman');
%             set(gca, 'YTick', [-1 0 1]);
%             %set(gca, 'XTick', [-i:1:-1]);
%             set(gca, 'TickLength', [0.015 0.025]);
            RowNo = RowNo + 1;
        end
    end

    FinalFRAxis = [min(FinalAxis(:,1)) max(FinalAxis(:,2)) 0 max(FinalAxis(:,4))];
    for i = 1:length(FRAxis),
        if (FRAxis(i) > 0)
            axes(FRAxis(i));
            axis(FinalFRAxis);
        end
    end
    
    CollapsedRasterFig = figure;
    set(gcf, 'Color', 'w');

    CollapsedINFRFig = figure;
    set(gcf, 'Color', 'w');
    
    CollapsedGapFRFig = figure;
    set(gcf, 'Color', 'w');
        
    for i = min(AllINPosition(:,end)):1:max(AllINPosition(:,end)),
        TempIndices = find(AllINPosition(:,end) == i);
        if (length(TempIndices) >= MinTrials)
            break;
        end
    end
    MinPosition = i;
    MinTrials = length(TempIndices);
    
    INNeuralAnalysisResults.ValidIndices = [];
    for i = MinPosition:1:max(AllINPosition(:,end)),
        INNeuralAnalysisResults.ValidIndices = [INNeuralAnalysisResults.ValidIndices; find(AllINPosition(:,end) == i)];
    end
    
    CollapsedINFR = [];
    CollapsedGapFR = [];
    
    PlotIndex = 0;
    for i = MinPosition:1:max(AllINPosition(:,end)),
        PlotIndex = PlotIndex + 1;
        PlotHGap = 0.03;
        PlotWGap = 0.075;

        PlotHt = 0.75;

        PlotWidth = 0.96/abs(MinPosition) - PlotWGap;

        TempRaster = [];
        TempIndices = find(AllINPosition(:,end) == i);
        rng('default')
        TempIndices = TempIndices(randperm(length(TempIndices)));
        
        for j = 1:length(TempIndices),
            TempRaster = [TempRaster; [AllINPreMotorSpikeTrain{TempIndices(j)} ones(size(AllINPreMotorSpikeTrain{TempIndices(j)}))*j]];
            TempRaster = [TempRaster; [AllGapPreMotorSpikeTrain{TempIndices(j)} ones(size(AllGapPreMotorSpikeTrain{TempIndices(j)}))*j]];
        end
                
        figure(CollapsedRasterFig);
        
        subplot('Position', [(0.04 + (abs(MinPosition) - abs(i))*(PlotWidth + PlotWGap)) 0.22 PlotWidth PlotHt]);
        hold on;
        
        if (~isempty(TempRaster))
            PlotRaster(TempRaster, 'k', 0.25, 0, MinTrials);
        end
        % plot([0 0], [-0.5 MinTrials+0.5], 'k');
        fill([(zeros(1,MinTrials) - PreMotorOffset) flipud(ones(1,MinTrials)*median(INDur) - PreMotorOffset)], [linspace(0.5, MinTrials+0.5, MinTrials) linspace(MinTrials+0.5, 0.5, MinTrials)], 'k', 'EdgeColor', 'none', 'FaceColor', [1 0.75 0.75], 'FaceAlpha', 0.5);
        fill([(ones(1,MinTrials)*median(INDur) - PreMotorOffset) flipud(ones(1,MinTrials)*median(INDur) - PreMotorOffset + LastGapDur)], [linspace(0.5, MinTrials+0.5, MinTrials) linspace(MinTrials+0.5, 0.5, MinTrials)], 'k', 'EdgeColor', 'none', 'FaceColor', [0.75 0.75 1], 'FaceAlpha', 0.5);
        if (strfind(PlotPSTH, 'on'))
            plot(INPSTEdges(:,1:end-1) + BinSize/2, mean(INPST{i}{j}(:,1:end-1))/max([max(AllINPST(:,1:end-1))]) * (MinTrials), 'Color', 'r', 'LineWidth', 1.5);
            plot(GapPSTEdges(:,1:end-1) + BinSize/2, mean(GapPST{i}{j}(:,1:end-1))/max([max(AllINPST(:,1:end-1))]) * (MinTrials), 'Color', 'b', 'LineWidth', 1.5);
        end
        set(gca, 'Layer', 'top');
        axis([-PreMotorOffset (median(INDur) - PreMotorOffset + LastGapDur) 0 MinTrials+0.5]);
        axis([-PreMotorOffset (median(INDur) - PreMotorOffset + LastGapDur) -1 MinTrials+0.5]);
        
        temp = axis;
        set(gca, 'XTick', [temp(1) 0], 'XTickLabel', [temp(1)*1000 0]);
        set(gca, 'Box', 'off');
        set(gca, 'FontSize', 8, 'FontName', 'Times New Roman');
        % set(gca, 'TickDir', 'out');
        set(gca, 'YTick', []);
        set(gca, 'YAxisLocation', 'right');
        set(gca, 'YColor', 'w');
        if (strfind(PlotPSTH, 'on'))
            if (j == 1)
                plot([-0.05 -0.05], [(MinTrials+0.5) ((MinTrials+0.5) - (100/max([max(AllINPST(:,1:end-1))]) * MinTrials))], 'r');
            end
        end
        CollapsedINFR = [CollapsedINFR; [ones(length(AllINPreMotorFR(TempIndices)),1)*i AllINPreMotorFR(TempIndices)]];
        CollapsedGapFR = [CollapsedGapFR; [ones(length(AllINGapPreMotorFR(TempIndices)),1)*i AllINGapPreMotorFR(TempIndices)]];
    end
    
    % Now to plot the firing rates
    figure(CollapsedINFRFig);
    hold on;
    set(gcf, 'Position', [920 277 150 420]);
    subplot('Position', [0.22 0.22 0.76 PlotHt]);
    hold on;
    for i = min(CollapsedINFR(:,1)):1:max(CollapsedINFR(:,1)),
        Indices = find(CollapsedINFR(:,1) == i);
        errorbar(mean(CollapsedINFR(Indices,1)), mean(CollapsedINFR(Indices,2)), std(CollapsedINFR(Indices,2)), 'ro', 'MarkerSize', 5); % firing rates using mean spike count in the window premotor to the IN
    end
    plot(unique(CollapsedINFR(:,1)), polyval(polyfit(CollapsedINFR(:,1), CollapsedINFR(:,2), 1), unique(CollapsedINFR(:,1))), 'r');
    [Rsq, F] = CalculateGoodnessofLinearFit(polyfit(CollapsedINFR(:,1), CollapsedINFR(:,2), 1), CollapsedINFR(:,1), CollapsedINFR(:,2));
    [r, p] = corrcoef(CollapsedINFR(:,1), CollapsedINFR(:,2));
    disp(['Rsq is ', num2str(Rsq), ' and corr. coeff ^2 is ', num2str(r(1,2).^2), ' and the p-value is ', num2str(p(1,2)), ' for INs']);
    axis tight;
    temp = axis;
    axis([(MinPosition-0.5) -0.75 0 temp(4)*1.1]);
    set(gca, 'Box', 'off');
    set(gca, 'XTick', []);
    set(gca, 'FontSize', 8, 'FontName', 'Times New Roman');
    %set(gca, 'YTickLabel', []);
    %set(gca, 'XTick', [-i:1:-1]);
    set(gca, 'TickLength', [0.015 0.025]);
    
    figure(CollapsedGapFRFig);
    hold on;
    set(gcf, 'Position', [920 277 150 420]);
    subplot('Position', [0.22 0.22 0.76 PlotHt]);
    hold on;
    MeanCollapsedGapFR = [];
    for i = min(CollapsedGapFR(:,1)):1:max(CollapsedGapFR(:,1)),
        Indices = find(CollapsedGapFR(:,1) == i);
        errorbar(mean(CollapsedGapFR(Indices,1)), mean(CollapsedGapFR(Indices,2)), std(CollapsedGapFR(Indices,2)), 'bo'); % firing rates using mean spike count in the window premotor to the IN
    end
    plot(unique(CollapsedGapFR(:,1)), polyval(polyfit(CollapsedGapFR(:,1), CollapsedGapFR(:,2), 1), unique(CollapsedGapFR(:,1))), 'b');
    [Rsq, F] = CalculateGoodnessofLinearFit(polyfit(CollapsedGapFR(:,1), CollapsedGapFR(:,2), 1), CollapsedGapFR(:,1), CollapsedGapFR(:,2));
    [r, p] = corrcoef(CollapsedGapFR(:,1), CollapsedGapFR(:,2));
    disp(['Rsq is ', num2str(Rsq), ' and corr. coeff ^2 is ', num2str(r(1,2).^2), ' and the p-value is ', num2str(p(1,2)), ' for Gaps']);
    axis tight;
    temp = axis;
    axis([(MinPosition-0.5) -0.75 0 temp(4)*1.1]);
    set(gca, 'Box', 'off');
    set(gca, 'XTick', []);
    set(gca, 'FontSize', 8, 'FontName', 'Times New Roman');
    %set(gca, 'YTickLabel', []);
    %set(gca, 'XTick', [-i:1:-1]);
    set(gca, 'TickLength', [0.015 0.025]);
    
    disp(['Firing rate scale bar is 100 Hz']);
end

disp('Finished Analysis');