function [INNeuralAnalysisResults] = IntroNoteNeuralAnalysisINSequence(Neural_INR, PlotOption, MinTrialNumber)

% Using a system for AllINPosition that keeps track of whether the intro
% note was the first, last or a middle intro note. The way I do this is by
% having 3 boolean flags for first, middle and last intro note. For
% instance if there was only one intro note, it would have the flags 1 0 1
% to indicate that is the first, it is also the last. 

% It also has 4 more flags - one for the total number of INs in that
% sequence, the second one gives the trial # for the particular sequence of
% INs and the third gives the position of the IN within that sequence and
% the fourth flag is the position of the intro note with reference to the
% last - the last being -1. This can also be obtained by doing (Column 6 -
% Column 4 - 1). With all of this it should be easy to reconstruct the 
% position and trial no of each IN.

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
AllINRaster = [];

MinNumber = 2;
BinSize = 0.01;
Edges = -0.05:BinSize:0.06;
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
                    AllINFeats(AllINIndex,:) = Neural_INR.BoutDetails(j).Feats(Neural_INR.INs{j}(k),1:4);
                    AllINFR(AllINIndex,:) = length(SpikeIndices)/(Edges(end) - Edges(1));
                    AllINIFR(AllINIndex,:) = interp1(Neural_INR.BoutDetails(j).IFR(1, IFRStartIndex:IFREndIndex), Neural_INR.BoutDetails(j).IFR(2, IFRStartIndex:IFREndIndex), StartTimes(k) + IFREdges);
                    AllINFR_IFR(AllINIndex,:) = mean(AllINIFR(AllINIndex,1:end-1));
                    if (~isempty(SpikeIndices))
                        AllINSpikeTrain{AllINIndex} = (Neural_INR.BoutDetails(j).SpikeTimes(SpikeIndices) - StartTimes(k));
                        AllINPST(AllINIndex,:) = histc(AllINSpikeTrain{AllINIndex}, Edges)/BinSize;
                        AllINRaster = [AllINRaster; [AllINSpikeTrain{AllINIndex} ones(size(AllINSpikeTrain{AllINIndex}))*AllINIndex]];
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
                    AllINFeats(AllINIndex,:) = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).Feats(Neural_INR.WithinBoutINs{j}(k),1:4);
                    AllINFR(AllINIndex,:) = length(SpikeIndices)/(Edges(end) - Edges(1));
                    AllINIFR(AllINIndex,:) = interp1(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).IFR(1, IFRStartIndex:IFREndIndex), Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).IFR(2, IFRStartIndex:IFREndIndex), StartTimes(k) + IFREdges);
                    AllINFR_IFR(AllINIndex,:) = mean(AllINIFR(AllINIndex,1:end-1));
                    if (~isempty(SpikeIndices))
                        AllINSpikeTrain{AllINIndex} = (Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(j)).SpikeTimes(SpikeIndices) - StartTimes(k));
                        AllINPST(AllINIndex,:) = histc(AllINSpikeTrain{AllINIndex}, Edges)/BinSize;
                        AllINRaster = [AllINRaster; [AllINSpikeTrain{AllINIndex} ones(size(AllINSpikeTrain{AllINIndex}))*AllINIndex]];
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

figure;
plot(IFREdges(1:end-1), mean(AllINIFR(:,1:end-1)), 'k');
hold on;
Threshold = mean(mean(AllINIFR(:,1:end-1))) + std(mean(AllINIFR(:,1:end-1)));
plot([IFREdges(1) IFREdges(end)], [Threshold Threshold], 'k--');
MeanIFR = mean(AllINIFR(:,1:end-1));
UpperCrossingIndex = find(MeanIFR >= Threshold, 1, 'first');
if (isempty(UpperCrossingIndex))
    UpperCrossingIndex = 1;
end
LowerCrossingIndex = UpperCrossingIndex + find(MeanIFR(UpperCrossingIndex+1:end) <= Threshold, 1, 'first');
if (isempty(LowerCrossingIndex))
    LowerCrossingIndex = length(MeanIFR);
end
%EventWindow = [(IFREdges(round(mean([LowerCrossingIndex UpperCrossingIndex]))) - 0.015) (IFREdges(round(mean([LowerCrossingIndex UpperCrossingIndex]))) + 0.015)];
EventWindow = [IFREdges(UpperCrossingIndex)-0.01 IFREdges(LowerCrossingIndex)+0.01];
plot([EventWindow(1) EventWindow(1)], [0 max(MeanIFR)], 'r');
plot([EventWindow(2) EventWindow(2)], [0 max(MeanIFR)], 'r');

figure;
if (~isempty(AllINRaster))
    PlotRaster(AllINRaster, 'k');
end
hold on;
plot([EventWindow(1) EventWindow(1)], [0 AllINIndex], 'r');
plot([EventWindow(2) EventWindow(2)], [0 AllINIndex], 'r');

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
if (MinTrials > 10)
    MinTrials = 10;
end

Temp_INIFR = INIFR;

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

INNeuralAnalysisResults.AllINSpikeTrain = AllINSpikeTrain;
INNeuralAnalysisResults.AllINPosition = AllINPosition;
INNeuralAnalysisResults.AllINIFR = AllINIFR;
INNeuralAnalysisResults.AllINFR = AllINFR;
INNeuralAnalysisResults.AllINFR_IFR = AllINFR_IFR;
INNeuralAnalysisResults.AllINPST = AllINPST;

if (exist('Corr', 'var'))
    INNeuralAnalysisResults.Corr = Corr;
    INNeuralAnalysisResults.StdCorr = StdCorr;
else
    INNeuralAnalysisResults.Corr = [];
    INNeuralAnalysisResults.StdCorr = [];
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

        PlotHt = 0.875/length(Indices) - PlotHGap;

        PlotWidth = 0.975/max(Indices) - PlotWGap;

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
                hold on;
                plot(Edges(1:end-1), mean(INPST{i}{j}(:,1:end-1))/max([mean(AllINPST(:,1:end-1))]) * MinTrials, 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
                if (~isempty(INRaster{i}{j}))
                    PlotRaster(INRaster{i}{j}, 'k', 0.25, 0, MinTrials);
                end
                fill([zeros(1,MinTrials) flipud(SyllOffsetRaster{i}{j}(1:MinTrials,1))'], [linspace(0.5, MinTrials+0.5, MinTrials) linspace(MinTrials+0.5, 0.5, MinTrials)], 'k', 'EdgeColor', 'none', 'FaceColor', [0.65 0.65 0.65], 'FaceAlpha', 0.5);
                plot(PrevSyllOffsetRaster{i}{j}(1:MinTrials,1), PrevSyllOffsetRaster{i}{j}(1:MinTrials,2), 'ro', 'MarkerSize', 4, 'MarkerFace', 'r');
                %plot(SyllOffsetRaster{i}{j}(1:MinTrials,1), SyllOffsetRaster{i}{j}(1:MinTrials,2), 'ro', 'MarkerSize', 2);
                plot(NextSyllOnsetRaster{i}{j}(1:MinTrials,1), NextSyllOnsetRaster{i}{j}(1:MinTrials,2), 'ro', 'MarkerSize', 4, 'MarkerFace', 'r');
                axis([Edges(1) Edges(end) 0 MinTrials+0.5]);
                %plot([0 0], [0.5 MinTrials+0.5], 'g--', 'LineWidth', 2);
                set(gca, 'YColor', 'w');

                if (j ~= 1)
                    set(gca, 'YTick', []);
                end

                if (i ~= Indices(end))
                    set(gca, 'XTick', []);
                    set(gca, 'XColor', 'w');
                else
                    if (j ~= round(max(Indices)/2))
                        set(gca, 'XTickLabel', []);
                    else
                        xlabel('Time (ms)', 'FontSize', 12, 'FontName', 'Arial');
                        set(gca, 'FontSize', 10, 'FontName', 'Arial');
                        set(gca, 'XTick', [-0.05 0 0.05], 'XTickLabel', [-50 0 50]);
                        set(gca, 'FontSize', 8, 'FontName', 'Arial');
                    end
                    set(gca, 'TickDir', 'out');
                    set(gca, 'TickLength', [0.075 0.25]);
                end
            end

            % Now to plot the firing rates
            figure(FRFig);
            set(gcf, 'Position', [920 277 150 420]);
            subplot('Position', [0.3 (1 - (RowNo+1)*(PlotHt + PlotHGap)) 0.65 PlotHt]);
            %errorbar([-i:1:-1], mean(INFR_IFR{i}(:,1:end-1)), zeros(size([-i:1:-1])), std(INFR_IFR{i}(:,1:end-1)), 'ks-'); % firing rates using mean IFR
            errorbar([-i:1:-1], mean(INFR{i}(:,1:end-1)), zeros(size([-i:1:-1])), std(INFR{i}(:,1:end-1)), 'ks-'); % firing rates using mean spike count
            axis tight;
            temp = axis;
            axis([(-max(Indices)-0.25) -0.75 0 temp(4)*1.1]);
            set(gca, 'Box', 'off');
            set(gca, 'FontSize', 8, 'FontName', 'Arial');
            if (i ~= max(Indices))
                if (i == Indices(round(length(Indices)/2)))
                    ylabel('Firing Rate (Hz)', 'FontSize', 10, 'FontName', 'Arial');
                end
                set(gca, 'XTickLabel', []);
            else
                set(gca, 'XTick', [-i:1:-1]);
                xlabel('Position of IN', 'FontSize', 12, 'FontName', 'Arial');
            end

            % Now to plot the correlations
            figure(CorrFig);
            set(gcf, 'Position', [1070 277 150 420]);
            subplot('Position', [0.3 (1 - (RowNo+1)*(PlotHt + PlotHGap)) 0.65 PlotHt]);
            errorbar([-i:1:-1], Corr{i}(1:end-1), zeros(size([-i:1:-1])), StdCorr{i}(1:end-1), 'ks-');
            hold on;
            errorbar([-i:1:-1], LastINCorr{i}(1:end-1), zeros(size([-i:1:-1])), LastINStdCorr{i}(1:end-1), 'rd-');
            axis tight;
            temp = axis;
            if (temp(4) > 1)
                axis([(-max(Indices)-0.25) -0.75 -1 temp(4)*1.1]);
            else
                axis([(-max(Indices)-0.25) -0.75 -1 1]);
            end
            set(gca, 'Box', 'off');
            set(gca, 'FontSize', 8, 'FontName', 'Arial');
            if (i ~= max(Indices))
                if (i == Indices(round(length(Indices)/2)))
                    ylabel('Pairwise correlation co-efficient', 'FontSize', 10, 'FontName', 'Arial');
                end
                set(gca, 'XTickLabel', []);
            else
                set(gca, 'XTick', [-i:1:-1]);
                xlabel('Position of IN', 'FontSize', 12, 'FontName', 'Arial');
            end
            RowNo = RowNo + 1;
        end
    end
    [coeff, score, latent] = princomp(AllINIFR(:,1:end-1));
    figure;
    hold on;
    plot(score(find(AllINPosition(:,1) == 1),1), score(find(AllINPosition(:,1) == 1),2), 'r+');
    plot(score(find(AllINPosition(:,3) == 1),1), score(find(AllINPosition(:,3) == 1),2), 'b+');
    plot(score(find(AllINPosition(:,2) == 1),1), score(find(AllINPosition(:,2) == 1),2), 'g+');
end

disp('Finished Analysis');