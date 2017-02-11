function [INNeuralAnalysisResults] = IntroNoteNeuralAnalysisCorrToFirstSyll(Neural_INR, BinSize)

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

PreMotorLag = 0.045;

Width = 0.01;
GaussianLen = 4;
IFRFs = 2000;
XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (IFRFs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * IFRFs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * IFRFs) * (Width * IFRFs)));


Edges = -1.5:BinSize:0.2;
INEdges = -0.1:BinSize:0.1;

IFREdges = INEdges(1):1/IFRFs:INEdges(end);

PST = [];
Index = 0;
TrialIndex = 0;

INIndex = 0;
AllPosition = [];

AllINRaster = [];
AllINPST = [];
AllINSpikeTimes = [];

AllGapRaster = [];
AllGapPST = [];
AllGapSpikeTimes = [];

AllINGapRaster = [];
AllINGapPST = [];
AllINGapSpikeTimes = [];

% First calculate the median IN duration across all INs and the median
% duration of all gaps between successive INs
% Also treat INs and gaps as one entity and calculate the median duration
% of the combined IN and gap durations

AllINDur = [];
AllGapDur = [];
AllINGapDur = [];

for i = 1:length(Neural_INR.NoofINs),
    if (Neural_INR.NoofINs(i) < 1)
        continue;
    end
    AllINDur = [AllINDur; (Neural_INR.BoutDetails(i).offsets(Neural_INR.INs{i}) - Neural_INR.BoutDetails(i).onsets(Neural_INR.INs{i}))];
    AllGapDur = [AllGapDur; (Neural_INR.BoutDetails(i).onsets(Neural_INR.INs{i} + 1) - Neural_INR.BoutDetails(i).offsets(Neural_INR.INs{i}))];
    AllINGapDur = [AllINGapDur; (Neural_INR.BoutDetails(i).onsets(Neural_INR.INs{i} + 1) - Neural_INR.BoutDetails(i).onsets(Neural_INR.INs{i}))];
end

for i = 1:size(Neural_INR.WithinBoutNoofINs,1),
    if (Neural_INR.WithinBoutNoofINs(i,1) < 1)
        continue;
    end
    AllINDur = [AllINDur; (Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).offsets(Neural_INR.WithinBoutINs{i}) - Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(Neural_INR.WithinBoutINs{i}))];
    AllGapDur = [AllGapDur; (Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(Neural_INR.WithinBoutINs{i} + 1) - Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).offsets(Neural_INR.WithinBoutINs{i}))];
    AllINGapDur = [AllINGapDur; (Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(Neural_INR.WithinBoutINs{i} + 1) - Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(Neural_INR.WithinBoutINs{i}))];
end

MedianINDur = median(AllINDur);
MedianGapDur = median(AllGapDur);
MedianINGapDur = median(AllINGapDur);

% Now find the spikes associated with each IN, each Gap and the combined
% IN-Gap for all INs and warp it to the median and plot the rasters by
% position

MaxINs = max([max(Neural_INR.NoofINs) max(Neural_INR.WithinBoutNoofINs(:,1))]);

for i = 1:MaxINs,
    INRaster{i} = [];
    INPST{i} = [];
    GapRaster{i} = [];
    GapPST{i} = [];
    INGapRaster{i} = [];
    INGapPST{i} = [];
end

INEdges = 0:BinSize:MedianINDur;
GapEdges = 0:BinSize:MedianGapDur;
INGapEdges = 0:BinSize:MedianINGapDur;


for i = 1:length(Neural_INR.NoofINs),
    if (Neural_INR.NoofINs(i) < 1)
        continue;
    end

    TrialIndex = TrialIndex + 1;
    
    for j = 1:Neural_INR.NoofINs(i),
        INStartTime = Neural_INR.BoutDetails(i).onsets(Neural_INR.INs{i}(j)) - PreMotorLag;
        INEndTime = Neural_INR.BoutDetails(i).offsets(Neural_INR.INs{i}(j)) - PreMotorLag;
        NextSyllStartTime = Neural_INR.BoutDetails(i).onsets(Neural_INR.INs{i}(j)+1) - PreMotorLag;

        INSpikeIndices = find((Neural_INR.BoutDetails(i).SpikeTimes >  (INStartTime)) & (Neural_INR.BoutDetails(i).SpikeTimes <=  (INEndTime)));
        INSpikeTimes{TrialIndex}{j} = Neural_INR.BoutDetails(i).SpikeTimes(INSpikeIndices) - INStartTime;
        INDur{TrialIndex}{j} = INEndTime - INStartTime;
        INWarpedSpikeTimes{TrialIndex}{j} = INSpikeTimes{TrialIndex}{j} * MedianINDur/INDur{TrialIndex}{j};
        
        INRaster{Neural_INR.NoofINs(i) - j + 1} = [INRaster{Neural_INR.NoofINs(i) - j + 1}; [INWarpedSpikeTimes{TrialIndex}{j} ones(size(INWarpedSpikeTimes{TrialIndex}{j}))*TrialIndex]];
        TempPST = histc(INWarpedSpikeTimes{TrialIndex}{j}, INEdges);
        if (size(TempPST,1) > size(TempPST,2))
            TempPST = TempPST';
        end
        INPST{Neural_INR.NoofINs(i) - j + 1} = [INPST{Neural_INR.NoofINs(i) - j + 1}; TempPST];
    
        GapSpikeIndices = find((Neural_INR.BoutDetails(i).SpikeTimes >  (INEndTime)) & (Neural_INR.BoutDetails(i).SpikeTimes <=  (NextSyllStartTime)));
        GapSpikeTimes{TrialIndex}{j} = Neural_INR.BoutDetails(i).SpikeTimes(GapSpikeIndices) - INEndTime;
        GapDur{TrialIndex}{j} = NextSyllStartTime - INEndTime;
        GapWarpedSpikeTimes{TrialIndex}{j} = GapSpikeTimes{TrialIndex}{j} * MedianGapDur/GapDur{TrialIndex}{j};
        
        GapRaster{Neural_INR.NoofINs(i) - j + 1} = [GapRaster{Neural_INR.NoofINs(i) - j + 1}; [GapWarpedSpikeTimes{TrialIndex}{j} ones(size(GapWarpedSpikeTimes{TrialIndex}{j}))*TrialIndex]];
        TempPST = histc(GapWarpedSpikeTimes{TrialIndex}{j}, GapEdges);
        if (size(TempPST,1) > size(TempPST,2))
            TempPST = TempPST';
        end
        GapPST{Neural_INR.NoofINs(i) - j + 1} = [GapPST{Neural_INR.NoofINs(i) - j + 1}; TempPST];
        
        INGapSpikeIndices = find((Neural_INR.BoutDetails(i).SpikeTimes >  (INStartTime)) & (Neural_INR.BoutDetails(i).SpikeTimes <=  (NextSyllStartTime)));
        INGapSpikeTimes{TrialIndex}{j} = Neural_INR.BoutDetails(i).SpikeTimes(INGapSpikeIndices) - INStartTime;
        INGapDur{TrialIndex}{j} = NextSyllStartTime - INStartTime;
        INGapWarpedSpikeTimes{TrialIndex}{j} = INGapSpikeTimes{TrialIndex}{j} * MedianINGapDur/INGapDur{TrialIndex}{j};
        
        INGapRaster{Neural_INR.NoofINs(i) - j + 1} = [INGapRaster{Neural_INR.NoofINs(i) - j + 1}; [INGapWarpedSpikeTimes{TrialIndex}{j} ones(size(INGapWarpedSpikeTimes{TrialIndex}{j}))*TrialIndex]];
        TempPST = histc(INGapWarpedSpikeTimes{TrialIndex}{j}, INGapEdges);
        if (size(TempPST,1) > size(TempPST,2))
            TempPST = TempPST';
        end
        INGapPST{Neural_INR.NoofINs(i) - j + 1} = [INGapPST{Neural_INR.NoofINs(i) - j + 1}; TempPST];
    end
end

for i = 1:size(Neural_INR.WithinBoutNoofINs,1),
    if (Neural_INR.WithinBoutNoofINs(i,1) < 1)
        continue;
    end
    
    Index = Index + 1;
    MotifStartTime = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(Neural_INR.WithinBoutINs{i}(end) + 1);
    SpikeIndices = find((Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes >  (MotifStartTime + Edges(1))) & (Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes <=  (MotifStartTime + Edges(end))));
    SpikeTimes = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes(SpikeIndices) - MotifStartTime;
    PST(Index,:) = histc(SpikeTimes, Edges);
    for j = 1:Neural_INR.WithinBoutNoofINs(i,1),
        INIndex = INIndex + 1;
        INStartTime = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(Neural_INR.WithinBoutINs{i}(j));
        SpikeIndices = find((Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes >  (INStartTime + INEdges(1))) & (Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes <=  (INStartTime + INEdges(end))));
        SpikeTimes = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes(SpikeIndices) - INStartTime;
        AllINPST(INIndex,:) = histc(SpikeTimes, INEdges);
        AllINPosition(INIndex) = j - Neural_INR.WithinBoutNoofINs(i,1) - 1;
        AllINDur(INIndex) = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).offsets(Neural_INR.WithinBoutINs{i}(j)) - Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(Neural_INR.WithinBoutINs{i}(j));
        AllINNextSyllStartTime(INIndex) = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(Neural_INR.WithinBoutINs{i}(j)+1) - INStartTime;
        
        IFRStartIndex = find(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).IFR(1,:) <= (INStartTime + Edges(1)), 1, 'last');
        IFREndIndex = find(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).IFR(1,:) <= (INStartTime + Edges(end)), 1, 'last');
        AllINIFR(INIndex,:) = interp1(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).IFR(1, IFRStartIndex:IFREndIndex), Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).IFR(2, IFRStartIndex:IFREndIndex), INStartTime + IFREdges);
    end
    INIndex = INIndex + 1;
    INStartTime = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(Neural_INR.WithinBoutINs{i}(end) + 1);
    SpikeIndices = find((Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes >  (INStartTime + INEdges(1))) & (Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes <=  (INStartTime + INEdges(end))));
    SpikeTimes = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).SpikeTimes(SpikeIndices) - INStartTime;
    AllINPST(INIndex,:) = histc(SpikeTimes, INEdges);
    AllINPosition(INIndex) = 0;
    AllINDur(INIndex) = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).offsets(Neural_INR.WithinBoutINs{i}(end) + 1) - Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(Neural_INR.WithinBoutINs{i}(end) + 1);
    AllINNextSyllStartTime(INIndex) = Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).onsets(Neural_INR.WithinBoutINs{i}(end) + 2) - INStartTime;
    IFRStartIndex = find(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).IFR(1,:) <= (INStartTime + Edges(1)), 1, 'last');
    IFREndIndex = find(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).IFR(1,:) <= (INStartTime + Edges(end)), 1, 'last');
    AllINIFR(INIndex,:) = interp1(Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).IFR(1, IFRStartIndex:IFREndIndex), Neural_INR.BoutDetails(Neural_INR.WithinBoutINBoutIndices(i)).IFR(2, IFRStartIndex:IFREndIndex), INStartTime + IFREdges);
end

INNeuralAnalysisResults.PST = PST;
INNeuralAnalysisResults.AllINPST = AllINPST;
INNeuralAnalysisResults.AllINPosition = AllINPosition;
INNeuralAnalysisResults.AllINDur = AllINDur;
INNeuralAnalysisResults.AllINIFR = AllINIFR;
INNeuralAnalysisResults.AllINNextSyllStartTime = AllINNextSyllStartTime;

disp('Finished Analysis');