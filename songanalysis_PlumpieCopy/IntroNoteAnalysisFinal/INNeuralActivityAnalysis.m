function [NeuralAnalysisResults] = INNeuralActivityAnalysis(NeuralResults)

MinNumber = 2;
Fs = 32000;
PSTEdges = 0.005;
PreDur = 0.05; % in seconds
PostDur = 0.015; % in seconds
FF_WindowWidth = 0.03; % Fano factor window width in seconds

% First calculate instantaneous firing rate for each bout
for i = 1:length(NeuralResults.BoutDetails),
    Time{i} = NeuralResults.BoutDetails(i).BoutOnset:1/Fs:NeuralResults.BoutDetails(i).BoutOffset;
    [IFR{i}] = CalculateIFR(NeuralResults.BoutDetails(i).SpikeTimes - NeuralResults.BoutDetails(i).BoutOnset, Time{i});
end
% First calculate activity at the onset of the bout depending on whether it
% is an IN-bout or a Non-IN bout

BoutIndices = find(NeuralResults.NoofINs > 0);

INBoutOnsetRaster = [];
INBoutNoofINs = NeuralResults.NoofINs(BoutIndices);

MaxINs = max(NeuralResults.NoofINs),

for i = 1:MaxINs,
    FFPST_FirstIN{i} = [];
    FFPST_SecondLastIN{i} = [];
    FFPST_LastIN{i} = [];
    FFPST_Motif{i} = [];
    for j = 1:length(BoutIndices),
        if (NeuralResults.NoofINs(BoutIndices(j)) == i)
            TrialSpikeTimes = [];
            TrialSpikeTimes = (NeuralResults.BoutDetails(BoutIndices(j)).SpikeTimes - NeuralResults.BoutDetails(BoutIndices(j)).onsets(NeuralResults.INs{BoutIndices(j)}(1)));
            ColIndex = 0;
            TempFFPST = [];
            for k = -1.5:0.001:0.2,
                ColIndex = ColIndex + 1;
                TempFFPST(1,ColIndex) = length(find((TrialSpikeTimes > k) & (TrialSpikeTimes <= (k+FF_WindowWidth))));
            end
            FFPST_FirstIN{i} = [FFPST_FirstIN{i}; TempFFPST];
            
            TrialSpikeTimes = [];
            TrialSpikeTimes = (NeuralResults.BoutDetails(BoutIndices(j)).SpikeTimes - NeuralResults.BoutDetails(BoutIndices(j)).onsets(NeuralResults.INs{BoutIndices(j)}(end-1)));
            ColIndex = 0;
            TempFFPST = [];
            for k = -1.5:0.001:0.2,
                ColIndex = ColIndex + 1;
                TempFFPST(1,ColIndex) = length(find((TrialSpikeTimes > k) & (TrialSpikeTimes <= (k+FF_WindowWidth))));
            end
            FFPST_SecondLastIN{i} = [FFPST_SecondLastIN{i}; TempFFPST];
            
            TrialSpikeTimes = [];
            TrialSpikeTimes = (NeuralResults.BoutDetails(BoutIndices(j)).SpikeTimes - NeuralResults.BoutDetails(BoutIndices(j)).onsets(NeuralResults.INs{BoutIndices(j)}(end)));
            ColIndex = 0;
            TempFFPST = [];
            for k = -1.5:0.001:0.2,
                ColIndex = ColIndex + 1;
                TempFFPST(1,ColIndex) = length(find((TrialSpikeTimes > k) & (TrialSpikeTimes <= (k+FF_WindowWidth))));
            end
            FFPST_LastIN{i} = [FFPST_LastIN{i}; TempFFPST];
            
            TrialSpikeTimes = [];
            TrialSpikeTimes = (NeuralResults.BoutDetails(BoutIndices(j)).SpikeTimes - NeuralResults.BoutDetails(BoutIndices(j)).onsets(NeuralResults.INs{BoutIndices(j)}(1) + 1));
            ColIndex = 0;
            TempFFPST = [];
            for k = -1.5:0.001:0.2,
                ColIndex = ColIndex + 1;
                TempFFPST(1,ColIndex) = length(find((TrialSpikeTimes > k) & (TrialSpikeTimes <= (k+FF_WindowWidth))));
            end
            FFPST_Motif{i} = [FFPST_Motif{i}; TempFFPST];
        end
    end
end
for i = 1:length(BoutIndices),
    Indices = find(NeuralResults.BoutDetails(i).SpikeTimes < NeuralResults.BoutDetails(i).onsets(1));
    TrialSpikeTimes = [];
    TrialSpikeTimes = (NeuralResults.BoutDetails(i).SpikeTimes(Indices,1) - NeuralResults.BoutDetails(i).onsets(1));
    INBoutOnsetRaster = [INBoutOnsetRaster; [TrialSpikeTimes ones(length(Indices),1)*i]];
    
    TimeIndex = find(Time{i} <= NeuralResults.BoutDetails(i).onsets(1), 1, 'last');
    INBoutOnsetIFR(i,:) = IFR{i}(1:TimeIndex);
    TrialSpikeTimes = (NeuralResults.BoutDetails(i).SpikeTimes - NeuralResults.BoutDetails(i).onsets(1));
    ColIndex = 0;
    for j = -1.5:0.001:0.05,
        ColIndex = ColIndex + 1;
        FFPST(i,ColIndex) = length(find((TrialSpikeTimes > j) & (TrialSpikeTimes <= (j+0.05))));
    end
end

BoutIndices = find((NeuralResults.FirstSyll ~= 'i') & (NeuralResults.SyllBeforeMotifs == 'i'));

Non_INBoutOnsetRaster = [];
Non_INBoutNoofINs = NeuralResults.NoofINs(BoutIndices);
for i = 1:length(BoutIndices),
    Indices = find(NeuralResults.BoutDetails(i).SpikeTimes < NeuralResults.BoutDetails(i).onsets(1));
    TrialSpikeTimes = [];
    TrialSpikeTimes = (NeuralResults.BoutDetails(i).SpikeTimes(Indices,1) - NeuralResults.BoutDetails(i).onsets(1));
    Non_INBoutOnsetRaster = [Non_INBoutOnsetRaster; [TrialSpikeTimes ones(length(Indices),1)*i]];
    
    TimeIndex = find(Time{i} <= NeuralResults.BoutDetails(i).onsets(1), 1, 'last');
    Non_INBoutOnsetIFR(i,:) = IFR{i}(1:TimeIndex);
end

BoutIndex = 0;
FiringRate = [];
SpikeTiming = [];
for i = 1:length(NeuralResults.NoofINs),
    BoutIndex = BoutIndex + 1;
    for j = 1:length(NeuralResults.INs{i}),
        if (j == 1)
            INPos = 2;
        else
            if (j == length(NeuralResults.INs{i}))
                INPos = 0;
            else
                INPos = 1;
            end
        end
        INStartTime = NeuralResults.BoutDetails(i).onsets(NeuralResults.INs{i}(j));
        INEndTime = NeuralResults.BoutDetails(i).offsets(NeuralResults.INs{i}(j));
        
        TrialSpikeTimes = [];
        TrialSpikeTimes = find((NeuralResults.BoutDetails(i).SpikeTimes >= (INStartTime - PreDur)) & (NeuralResults.BoutDetails(i).SpikeTimes <= (INEndTime + PostDur)));
              
        IFRStartIndex = find(Time{i} <= (INStartTime - PreDur), 1, 'last');
        IFREndIndex = find(Time{i} <= (INEndTime + PostDur), 1, 'last');
        
        IN_IFR{BoutIndex}{j} = [(Time{i}(IFRStartIndex:IFREndIndex) - INStartTime); IFR{i}(IFRStartIndex:IFREndIndex)];
        FiringRate = [FiringRate; [mean(IN_IFR{BoutIndex}{j}(2,:)) (length(NeuralResults.INs{i}) - j + 1) INPos]];
        if (~isempty(TrialSpikeTimes))
            Raster{BoutIndex}{j} = NeuralResults.BoutDetails(i).SpikeTimes(TrialSpikeTimes) - INStartTime;
            SpikeTiming = [SpikeTiming; [Raster{BoutIndex}{j}(1) (length(NeuralResults.INs{i}) - j + 1) INPos]];
        end
    end
end

BoutIndices = find(NeuralResults.WithinBoutNoofINs(:,1) > 0);

for i = 1:length(BoutIndices),
    BoutIndex = BoutIndex + 1;
    for j = 1:length(NeuralResults.WithinBoutINs{BoutIndices(i)}),
        if (j == 1)
            INPos = 2;
        else
            if (j == length(NeuralResults.WithinBoutINs{BoutIndices(i)}))
                INPos = 0;
            else
                INPos = 1;
            end
        end
        INStartTime = NeuralResults.BoutDetails(NeuralResults.WithinBoutINBoutIndices(BoutIndices(i))).onsets(NeuralResults.WithinBoutINs{BoutIndices(i)}(j));
        INEndTime = NeuralResults.BoutDetails(NeuralResults.WithinBoutINBoutIndices(BoutIndices(i))).offsets(NeuralResults.WithinBoutINs{BoutIndices(i)}(j));
        
        TrialSpikeTimes = [];
        TrialSpikeTimes = find((NeuralResults.BoutDetails(NeuralResults.WithinBoutINBoutIndices(BoutIndices(i))).SpikeTimes >= (INStartTime - PreDur)) & (NeuralResults.BoutDetails(NeuralResults.WithinBoutINBoutIndices(BoutIndices(i))).SpikeTimes <= (INEndTime + PostDur)));
                
        IFRStartIndex = find(Time{NeuralResults.WithinBoutINBoutIndices(BoutIndices(i))} <= (INStartTime - PreDur), 1, 'last');
        IFREndIndex = find(Time{NeuralResults.WithinBoutINBoutIndices(BoutIndices(i))} <= (INEndTime + PostDur), 1, 'last');
        
        IN_IFR{BoutIndex}{j} = [(Time{NeuralResults.WithinBoutINBoutIndices(BoutIndices(i))}(IFRStartIndex:IFREndIndex) - INStartTime); IFR{NeuralResults.WithinBoutINBoutIndices(BoutIndices(i))}(IFRStartIndex:IFREndIndex)];
        FiringRate = [FiringRate; [mean(IN_IFR{BoutIndex}{j}(2,:)) (length(NeuralResults.WithinBoutINs{BoutIndices(i)}) - j + 1) INPos]];
        if (~isempty(TrialSpikeTimes))
            Raster{BoutIndex}{j} = NeuralResults.BoutDetails(NeuralResults.WithinBoutINBoutIndices(BoutIndices(i))).SpikeTimes(TrialSpikeTimes) - INStartTime;
            SpikeTiming = [SpikeTiming; [Raster{BoutIndex}{j}(1) (length(NeuralResults.WithinBoutINs{BoutIndices(i)}) - j + 1) INPos]];
        end
    end
end

disp('Finished neural analysis');
