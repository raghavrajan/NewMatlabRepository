function [UnDirBoutSpikingDetails, DirBoutSpikingDetails] = LSINA_CalcSpikeCountRasterPST(SURecordingDetails)

% This function will calculate the spontaneous activity, the pre-bout
% activity, the pst, the raster and the first syllable for all bouts

SpontActivityWindow = [-2000 -1500];
PreSongActivityWindow = [-600 -100];
PreMotorSongActivityWindow = [-100 0];

PSTBinSize = 0.010; % in sec
PSTEdges = -2000:PSTBinSize*1000:200; 

SmoothFR_TimeStep = 0.1; % in ms

% Gaussian window for smoothing
GaussianLen = 1;
Width = 0.025; % in sec
XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (1/PSTBinSize)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * 1/PSTBinSize) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * 1/PSTBinSize) * (Width * 1/PSTBinSize)));
SmoothGaussWin = (1/((Width * 1/(SmoothFR_TimeStep/1000)) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * 1/(SmoothFR_TimeStep/1000)) * (Width * 1/(SmoothFR_TimeStep/1000))));
InterBoutInterval = SURecordingDetails.Interboutinterval;

% First load up all the data using a file that does the pre processing

disp('Getting song, non-song rasters, psts, etc. ...');
% If this is continuous data, then put all spikes together in one long
% array.
if (SURecordingDetails.Continuousdata == 1)

    % Put all syllable labels together in one long array
    [AllLabels, AllOnsets, AllOffsets] = CombineContinuousDataNoteFiles(SURecordingDetails.DataDirectory, SURecordingDetails.SongFileNames, fullfile(SURecordingDetails.DataDirectory, 'ASSLNoteFiles'), SURecordingDetails.FileType);
     % Now get rid of female in, out events and female calls
    FemaleINOUTEvents = regexp(AllLabels, ['[', SURecordingDetails.Femalecalllabel, SURecordingDetails.FemaleinLabel, SURecordingDetails.FemaleoutLabel, SURecordingDetails.Closedoorlabel, ']']);
    if (~isempty(FemaleINOUTEvents))
        AllLabels(FemaleINOUTEvents) = [];
        AllOnsets(FemaleINOUTEvents) = [];
        AllOffsets(FemaleINOUTEvents) = [];
    end

    ClusterSpikeData = [];
    AllClusterSpikeData = []; % includes optional cluster times as well
    TotalFileLen = 0;

    % First get all spike times into one long array
    for j = 1:length(SURecordingDetails.ClusterSpikeTimes),
        if (~isempty(SURecordingDetails.ClusterSpikeTimes{j}))
            ClusterSpikeData = [ClusterSpikeData; (TotalFileLen + SURecordingDetails.ClusterSpikeTimes{j}(:))];
        end
        if (~isempty(SURecordingDetails.OptionalClusterSpikeTimes{j}))
            AllClusterSpikeData = [AllClusterSpikeData; (TotalFileLen + SURecordingDetails.OptionalClusterSpikeTimes{j}(:))];
        end
        TotalFileLen = TotalFileLen + SURecordingDetails.FileLen(j);
    end
end

% First find all undir bouts that have more than the inter-bout
% interval before the start and at the end
% 25/02/2017 - realised a few days ago that the 9th column of Bouts has to
% be > 1 for the bout to considered as having enough post-data - have to
% redo some of the analysis and figures because of this.

UnDirBouts = find((SURecordingDetails.BoutDirUnDir == 0) & (SURecordingDetails.Bouts(:,8) > 0)' & (SURecordingDetails.Bouts(:,9) > 1)');

if (isempty(UnDirBouts))
   UnDirBoutSpikingDetails = [];
end

% Now for each of the neurons, take the bouts with enough data at the
% beginning and find all the spikes before bout initiation in two different
% windows - -2 to -1.5s (spont. activity window) and -0.6 to -0.1s
% (pre-song window)
% In addition get the spike train, the raster and the pst for that bout.
% Keep track of bout # also

if (SURecordingDetails.Continuousdata == 1)
    Index = 1;

    for j = UnDirBouts(:)',
        
        UnDirBoutSpikingDetails.BoutIndex(Index) = j;
        if (SURecordingDetails.Bouts(j,3) == 1)
            BoutOnsetTime = SURecordingDetails.Bouts(j,5);
        else
            BoutOnsetTime = sum(SURecordingDetails.FileLen(1:(SURecordingDetails.Bouts(j,3) - 1))) + SURecordingDetails.Bouts(j,5);
        end
        if (SURecordingDetails.Bouts(j,4) == 1)
            BoutOffsetTime = SURecordingDetails.Bouts(j,6);
        else
            BoutOffsetTime = sum(SURecordingDetails.FileLen(1:(SURecordingDetails.Bouts(j,4) - 1))) + SURecordingDetails.Bouts(j,6);
        end
        
        UnDirBoutSpikingDetails.SmoothBoutSpikeTrain{Index} = zeros(length(-2000:SmoothFR_TimeStep:(BoutOffsetTime - BoutOnsetTime + 2000)), 1);
        
        UnDirBoutSpikingDetails.SpontWindowSpikeCount(Index) = length(find((ClusterSpikeData >= (BoutOnsetTime + SpontActivityWindow(1))) & (ClusterSpikeData < (BoutOnsetTime + SpontActivityWindow(2)))))/((SpontActivityWindow(2) - SpontActivityWindow(1))/1000);
        UnDirBoutSpikingDetails.PreBoutWindowSpikeCount(Index) = length(find((ClusterSpikeData >= (BoutOnsetTime + PreSongActivityWindow(1))) & (ClusterSpikeData < (BoutOnsetTime + PreSongActivityWindow(2)))))/((SpontActivityWindow(2) - SpontActivityWindow(1))/1000);
                
        UnDirBoutSpikingDetails.PreMotorBoutWindowSpikeCount(Index) = length(find((ClusterSpikeData >= (BoutOnsetTime + PreMotorSongActivityWindow(1))) & (ClusterSpikeData < (BoutOnsetTime + PreMotorSongActivityWindow(2)))))/((PreMotorSongActivityWindow(2) - PreMotorSongActivityWindow(1))/1000);
        UnDirBoutSpikingDetails.BoutOnsetSpikeTrain{Index} = ClusterSpikeData(find((ClusterSpikeData >= (BoutOnsetTime -2000)) & (ClusterSpikeData < (BoutOnsetTime + 200)))) - BoutOnsetTime;
        UnDirBoutSpikingDetails.BoutSpikeTrain{Index} = ClusterSpikeData(find((ClusterSpikeData >= (BoutOnsetTime - InterBoutInterval)) & (ClusterSpikeData <= (BoutOffsetTime + InterBoutInterval)))) - BoutOnsetTime;
        if (~isempty(UnDirBoutSpikingDetails.BoutSpikeTrain{Index}))
            for Spikes = 1:length(UnDirBoutSpikingDetails.BoutSpikeTrain{Index}),
                SpikeLocationIndex = round(UnDirBoutSpikingDetails.BoutSpikeTrain{Index}(Spikes) * 1 / SmoothFR_TimeStep);
                if (SpikeLocationIndex < 1)
                    SpikeLocationIndex = 1;
                else
                    if (SpikeLocationIndex > length(UnDirBoutSpikingDetails.SmoothBoutSpikeTrain{Index}))
                        SpikeLocationIndex = length(UnDirBoutSpikingDetails.SmoothBoutSpikeTrain{Index});
                    end
                end
                UnDirBoutSpikingDetails.SmoothBoutSpikeTrain{Index}(SpikeLocationIndex) = 1;
            end
            UnDirBoutSpikingDetails.SmoothBoutSpikeTrain{Index} = conv(UnDirBoutSpikingDetails.SmoothBoutSpikeTrain{Index}, SmoothGaussWin, 'same');
        end
    
                        
        
        if (~isempty(UnDirBoutSpikingDetails.BoutOnsetSpikeTrain{Index}))
            UnDirBoutSpikingDetails.SpikeRaster{Index} = [UnDirBoutSpikingDetails.BoutOnsetSpikeTrain{Index}(:) ones(size(UnDirBoutSpikingDetails.BoutOnsetSpikeTrain{Index}(:)))];
        else
            UnDirBoutSpikingDetails.SpikeRaster{Index} = [];
        end
        
        UnDirBoutSpikingDetails.BoutLabels{Index} = AllLabels(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        UnDirBoutSpikingDetails.BoutOnsets{Index} = AllOnsets(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        UnDirBoutSpikingDetails.BoutOffsets{Index} = AllOffsets(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        
        % Also find spikes in the 500ms before the first 'i' in the bout
        Flag = 0;
        for SyllIndex = 1:length(UnDirBoutSpikingDetails.BoutLabels{Index})
            if (~isempty(find(SURecordingDetails.INLabels == UnDirBoutSpikingDetails.BoutLabels{Index}(SyllIndex))))
                UnDirBoutSpikingDetails.FirstINIndex(Index) = SyllIndex;
                UnDirBoutSpikingDetails.FirstINOnsetTime(Index) = UnDirBoutSpikingDetails.BoutOnsets{Index}(SyllIndex);
                UnDirBoutSpikingDetails.PreFirstINPreBoutWindowSpikeCount(Index) = length(find((ClusterSpikeData >= (UnDirBoutSpikingDetails.FirstINOnsetTime(Index) + PreSongActivityWindow(1))) & (ClusterSpikeData < (UnDirBoutSpikingDetails.FirstINOnsetTime(Index) + PreSongActivityWindow(2)))))/((SpontActivityWindow(2) - SpontActivityWindow(1))/1000);
                Flag = 1;
                break;
            end
        end
        
        if (Flag == 0)
            UnDirBoutSpikingDetails.FirstINIndex(Index) = NaN;
            UnDirBoutSpikingDetails.FirstINOnsetTime(Index) = NaN;
            UnDirBoutSpikingDetails.PreFirstINPreBoutWindowSpikeCount(Index) = NaN;
        end
        
        UnDirBoutSpikingDetails.BoutFirstSyllLabel(Index) = UnDirBoutSpikingDetails.BoutLabels{Index}(1);
        UnDirBoutSpikingDetails.BoutPST(Index,:) = conv(histc(ClusterSpikeData, PSTEdges + BoutOnsetTime)/PSTBinSize, GaussWin, 'same');
        UnDirBoutSpikingDetails.SongBoutOrNot(Index) = SURecordingDetails.Bouts(j,7);
        
        Index = Index + 1;
    end
else
    Index = 1;
    for j = UnDirBouts(:)',
        BoutOnsetTime = SURecordingDetails.Bouts(j,5);
        BoutOffsetTime = SURecordingDetails.Bouts(j,6);
 
        UnDirBoutSpikingDetails.BoutIndex(Index) = j;
        
        UnDirBoutSpikingDetails.SpontWindowSpikeCount(Index) = length(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutOnsetTime + SpontActivityWindow(1))) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} < (BoutOnsetTime + SpontActivityWindow(2)))))/((SpontActivityWindow(2) - SpontActivityWindow(1))/1000);
        UnDirBoutSpikingDetails.PreBoutWindowSpikeCount(Index) = length(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutOnsetTime + PreSongActivityWindow(1))) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} < (BoutOnsetTime + PreSongActivityWindow(2)))))/((PreSongActivityWindow(2) - PreSongActivityWindow(1))/1000);
        UnDirBoutSpikingDetails.PreMotorBoutWindowSpikeCount(Index) = length(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutOnsetTime + PreMotorSongActivityWindow(1))) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} < (BoutOnsetTime + PreMotorSongActivityWindow(2)))))/((PreMotorSongActivityWindow(2) - PreMotorSongActivityWindow(1))/1000);
        UnDirBoutSpikingDetails.BoutOnsetSpikeTrain{Index} = SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)}(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutOnsetTime - 2000)) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} < (BoutOnsetTime + 200)))) - BoutOnsetTime;
        UnDirBoutSpikingDetails.BoutSpikeTrain{Index} = SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)}(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutOnsetTime - InterBoutInterval)) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} <= (BoutOffsetTime + InterBoutInterval)))) - BoutOnsetTime;
        
        if (~isempty(UnDirBoutSpikingDetails.BoutOnsetSpikeTrain{Index}))
            UnDirBoutSpikingDetails.SpikeRaster{Index} = [UnDirBoutSpikingDetails.BoutOnsetSpikeTrain{Index}(:) ones(size(UnDirBoutSpikingDetails.BoutOnsetSpikeTrain{Index}(:)))];
        else
            UnDirBoutSpikingDetails.SpikeRaster{Index} = [];
        end
        
        UnDirBoutSpikingDetails.BoutLabels{Index} = SURecordingDetails.NoteInfo{SURecordingDetails.Bouts(j,3)}.labels(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        UnDirBoutSpikingDetails.BoutOnsets{Index} = SURecordingDetails.NoteInfo{SURecordingDetails.Bouts(j,3)}.onsets(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        UnDirBoutSpikingDetails.BoutOffsets{Index} = SURecordingDetails.NoteInfo{SURecordingDetails.Bouts(j,3)}.offsets(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        
        % Also find spikes in the 500ms before the first 'i' in the bout
        Flag = 0;
        for SyllIndex = 1: length(UnDirBoutSpikingDetails.BoutLabels{Index})
            if (~isempty(find(SURecordingDetails.INLabels == UnDirBoutSpikingDetails.BoutLabels{Index}(SyllIndex))))
                UnDirBoutSpikingDetails.FirstINIndex(Index) = SyllIndex;
                UnDirBoutSpikingDetails.FirstINOnsetTime(Index) = UnDirBoutSpikingDetails.BoutOnsets{Index}(SyllIndex);
                UnDirBoutSpikingDetails.PreFirstINPreBoutWindowSpikeCount(Index) = length(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (UnDirBoutSpikingDetails.FirstINOnsetTime(Index) + PreSongActivityWindow(1))) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} < (UnDirBoutSpikingDetails.FirstINOnsetTime(Index) + PreSongActivityWindow(2)))))/((SpontActivityWindow(2) - SpontActivityWindow(1))/1000);
                Flag = 1;
                break;
            end
        end
        
        if (Flag == 0)
            UnDirBoutSpikingDetails.FirstINIndex(Index) = NaN;
            UnDirBoutSpikingDetails.FirstINOnsetTime(Index) = NaN;
            UnDirBoutSpikingDetails.PreFirstINPreBoutWindowSpikeCount(Index) = NaN;
        end
        
        UnDirBoutSpikingDetails.BoutFirstSyllLabel(Index) = UnDirBoutSpikingDetails.BoutLabels{Index}(1);
        UnDirBoutSpikingDetails.BoutPST(Index,:) = conv(histc(SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)}, PSTEdges + BoutOnsetTime)/PSTBinSize, GaussWin, 'same');
        UnDirBoutSpikingDetails.SongBoutOrNot(Index) = SURecordingDetails.Bouts(j,7);
        
        Index = Index + 1;
    end
end

% Now do the same for directed song - only for directed song I will also
% consider all bouts that have enough time in the pre-bout window. I can
% always compare everything to spontaneous activity in undir song.

DirBouts = find((SURecordingDetails.BoutDirUnDir == 1));

if (isempty(DirBouts))
   DirBoutSpikingDetails = [];
end

% Now for each of the neurons, take the bouts with enough data at the
% beginning and find all the spikes before bout initiation in two different
% windows - -2 to -1.5s (spont. activity window) and -0.6 to -0.1s
% (pre-song window)
% In addition get the spike train, the raster and the pst for that bout.
% Keep track of bout # also

if (SURecordingDetails.Continuousdata == 1)
    Index = 1;

    for j = DirBouts(:)',
        
        DirBoutSpikingDetails.BoutIndex(Index) = j;
        if (SURecordingDetails.Bouts(j,3) == 1)
            BoutOnsetTime = SURecordingDetails.Bouts(j,5);
        else
            BoutOnsetTime = sum(SURecordingDetails.FileLen(1:(SURecordingDetails.Bouts(j,3) - 1))) + SURecordingDetails.Bouts(j,5);
        end
        if (SURecordingDetails.Bouts(j,4) == 1)
            BoutOffsetTime = SURecordingDetails.Bouts(j,6);
        else
            BoutOffsetTime = sum(SURecordingDetails.FileLen(1:(SURecordingDetails.Bouts(j,4) - 1))) + SURecordingDetails.Bouts(j,6);
        end
        
        if (BoutOnsetTime >= InterBoutInterval)
            DirBoutSpikingDetails.SpontWindowSpikeCount(Index) = length(find((ClusterSpikeData >= (BoutOnsetTime + SpontActivityWindow(1))) & (ClusterSpikeData < (BoutOnsetTime + SpontActivityWindow(2)))))/((SpontActivityWindow(2) - SpontActivityWindow(1))/1000);
            DirBoutSpikingDetails.BoutOnsetSpikeTrain{Index} = ClusterSpikeData(find((ClusterSpikeData >= (BoutOnsetTime -2000)) & (ClusterSpikeData < (BoutOnsetTime + 200)))) - BoutOnsetTime;
            DirBoutSpikingDetails.BoutPST(Index,:) = conv(histc(ClusterSpikeData, PSTEdges + BoutOnsetTime)/PSTBinSize, GaussWin, 'same');
            if (BoutOffsetTime <= (sum(SURecordingDetails.FileLen) - InterBoutInterval))
                DirBoutSpikingDetails.BoutSpikeTrain{Index} = ClusterSpikeData(find((ClusterSpikeData >= (BoutOnsetTime - InterBoutInterval)) & (ClusterSpikeData <= (BoutOffsetTime + InterBoutInterval)))) - BoutOnsetTime;
            else
                DirBoutSpikingDetails.BoutSpikeTrain{Index} = [];
            end
        else
            DirBoutSpikingDetails.BoutOnsetSpikeTrain{Index} = [];
            DirBoutSpikingDetails.BoutPST(Index,:) = ones(1,length(PSTEdges))*NaN;
            DirBoutSpikingDetails.SpontWindowSpikeCount(Index) = -1; % a negative value to indicate that there was no data
            DirBoutSpikingDetails.BoutSpikeTrain{Index} = [];
        end
        
        if (BoutOnsetTime >= PreSongActivityWindow(1))
            DirBoutSpikingDetails.PreBoutWindowSpikeCount(Index) = length(find((ClusterSpikeData >= (BoutOnsetTime + PreSongActivityWindow(1))) & (ClusterSpikeData < (BoutOnsetTime + PreSongActivityWindow(2)))))/((SpontActivityWindow(2) - SpontActivityWindow(1))/1000);
            DirBoutSpikingDetails.PreMotorBoutWindowSpikeCount(Index) = length(find((ClusterSpikeData >= (BoutOnsetTime + PreMotorSongActivityWindow(1))) & (ClusterSpikeData < (BoutOnsetTime + PreMotorSongActivityWindow(2)))))/((PreMotorSongActivityWindow(2) - PreMotorSongActivityWindow(1))/1000);
            DirBoutSpikingDetails.BoutOnsetSpikeTrain{Index} = ClusterSpikeData(find((ClusterSpikeData >= (BoutOnsetTime -2000)) & (ClusterSpikeData < (BoutOnsetTime + 200)))) - BoutOnsetTime;
            DirBoutSpikingDetails.BoutPST(Index,:) = conv(histc(ClusterSpikeData, PSTEdges + BoutOnsetTime)/PSTBinSize, GaussWin, 'same');
        else
            DirBoutSpikingDetails.PreBoutWindowSpikeCount(Index) = -1;
        end
        
        if (~isempty(DirBoutSpikingDetails.BoutOnsetSpikeTrain{Index}))
            DirBoutSpikingDetails.SpikeRaster{Index} = [DirBoutSpikingDetails.BoutOnsetSpikeTrain{Index}(:) ones(size(DirBoutSpikingDetails.BoutOnsetSpikeTrain{Index}(:)))];
        else
            DirBoutSpikingDetails.SpikeRaster{Index} = [];
        end
        
        DirBoutSpikingDetails.BoutLabels{Index} = AllLabels(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        DirBoutSpikingDetails.BoutOnsets{Index} = AllOnsets(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        DirBoutSpikingDetails.BoutOffsets{Index} = AllOffsets(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        DirBoutSpikingDetails.BoutFirstSyllLabel(Index) = DirBoutSpikingDetails.BoutLabels{Index}(1);
        
        DirBoutSpikingDetails.SongBoutOrNot(Index) = SURecordingDetails.Bouts(j,7);
        
        Index = Index + 1;
    end
else
    Index = 1;
    for j = DirBouts(:)',
        BoutOnsetTime = SURecordingDetails.Bouts(j,5);
        BoutOffsetTime = SURecordingDetails.Bouts(j,6);
 
        DirBoutSpikingDetails.BoutIndex(Index) = j;
        
        if (BoutOnsetTime >= InterBoutInterval)
            DirBoutSpikingDetails.SpontWindowSpikeCount(Index) = length(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutOnsetTime + SpontActivityWindow(1))) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} < (BoutOnsetTime + SpontActivityWindow(2)))))/((SpontActivityWindow(2) - SpontActivityWindow(1))/1000);
            if ((BoutOnsetTime + 200) <= SURecordingDetails.FileLen(SURecordingDetails.Bouts(j,3)))
                DirBoutSpikingDetails.BoutOnsetSpikeTrain{Index} = SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)}(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutOnsetTime - 2000)) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} < (BoutOnsetTime + 200)))) - BoutOnsetTime;
                DirBoutSpikingDetails.BoutPST(Index,:) = conv(histc(SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)}, PSTEdges + BoutOnsetTime)/PSTBinSize, GaussWin, 'same');
            else
                DirBoutSpikingDetails.BoutOnsetSpikeTrain{Index} = [];
                DirBoutSpikingDetails.BoutPST(Index,:) = ones(1, length(PSTEdges))*NaN;
            end
            if ((BoutOffsetTime + InterBoutInterval) <= SURecordingDetails.FileLen(SURecordingDetails.Bouts(j,3)))
                DirBoutSpikingDetails.BoutSpikeTrain{Index} = SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)}(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutOnsetTime - InterBoutInterval)) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} <= (BoutOffsetTime + InterBoutInterval)))) - BoutOnsetTime;
            else
                DirBoutSpikingDetails.BoutSpikeTrain{Index} = [];
            end
        else
            DirBoutSpikingDetails.BoutOnsetSpikeTrain{Index} = [];
            DirBoutSpikingDetails.BoutPST(Index,:) = ones(1,length(PSTEdges))*NaN;
            DirBoutSpikingDetails.SpontWindowSpikeCount(Index) = -1; % a negative value to indicate that there was no data
            DirBoutSpikingDetails.BoutSpikeTrain{Index} = [];
        end
        
        if (BoutOnsetTime >= PreSongActivityWindow(1))
            DirBoutSpikingDetails.PreBoutWindowSpikeCount(Index) = length(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutOnsetTime + PreSongActivityWindow(1))) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} < (BoutOnsetTime + PreSongActivityWindow(2)))))/((SpontActivityWindow(2) - SpontActivityWindow(1))/1000);
            DirBoutSpikingDetails.PreMotorBoutWindowSpikeCount(Index) = length(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutOnsetTime + PreMotorSongActivityWindow(1))) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} < (BoutOnsetTime + PreMotorSongActivityWindow(2)))))/((PreMotorSongActivityWindow(2) - PreMotorSongActivityWindow(1))/1000);
            DirBoutSpikingDetails.BoutOnsetSpikeTrain{Index} = SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)}(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutOnsetTime - 2000)) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} < (BoutOnsetTime + 200)))) - BoutOnsetTime;
            DirBoutSpikingDetails.BoutPST(Index,:) = conv(histc(SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)}, PSTEdges + BoutOnsetTime)/PSTBinSize, GaussWin, 'same');
            % In this case, bout onset spike train will have only the times
            % for the pre-song window. If spont. window spike count is -1,
            % then I know that the raster should be limited to the pre-song
            % window alone
        else
            DirBoutSpikingDetails.PreBoutWindowSpikeCount(Index) = -1;
        end
        
        if (~isempty(DirBoutSpikingDetails.BoutOnsetSpikeTrain{Index}))
            DirBoutSpikingDetails.SpikeRaster{Index} = [DirBoutSpikingDetails.BoutOnsetSpikeTrain{Index}(:) ones(size(DirBoutSpikingDetails.BoutOnsetSpikeTrain{Index}(:)))];
        else
            DirBoutSpikingDetails.SpikeRaster{Index} = [];
        end
        
        DirBoutSpikingDetails.BoutLabels{Index} = SURecordingDetails.NoteInfo{SURecordingDetails.Bouts(j,3)}.labels(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        DirBoutSpikingDetails.BoutOnsets{Index} = SURecordingDetails.NoteInfo{SURecordingDetails.Bouts(j,3)}.onsets(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        DirBoutSpikingDetails.BoutOffsets{Index} = SURecordingDetails.NoteInfo{SURecordingDetails.Bouts(j,3)}.offsets(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        
        DirBoutSpikingDetails.BoutFirstSyllLabel(Index) = DirBoutSpikingDetails.BoutLabels{Index}(1);
        
        DirBoutSpikingDetails.SongBoutOrNot(Index) = SURecordingDetails.Bouts(j,7);
        
        Index = Index + 1;
    end
end