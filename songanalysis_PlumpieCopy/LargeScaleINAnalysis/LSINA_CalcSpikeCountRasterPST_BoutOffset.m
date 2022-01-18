function [UnDirBoutSpikingDetails, DirBoutSpikingDetails] = LSINA_CalcSpikeCountRasterPST_BoutOffset(SURecordingDetails)

% This function will calculate the spontaneous activity, the pre-bout
% activity, the pst, the raster and the first syllable for all bouts

SpontActivityWindow = [1500 2000];
PostSongActivityWindow = [100 600];
PostMotorSongActivityWindow = [0 100];

PSTBinSize = 0.010; % in sec
PSTEdges = -200:PSTBinSize*1000:2000; 

% Gaussian window for smoothing
GaussianLen = 1;
Width = 0.025; % in sec
XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (1/PSTBinSize)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * 1/PSTBinSize) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * 1/PSTBinSize) * (Width * 1/PSTBinSize)));

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
% beginning and find all the spikes after bout end in two different
% windows - 1.5 to 2s (spont. activity window) and 0.1 to 0.6s
% (post-song window)
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
        
        UnDirBoutSpikingDetails.SpontWindowSpikeCount(Index) = length(find((ClusterSpikeData >= (BoutOffsetTime + SpontActivityWindow(1))) & (ClusterSpikeData < (BoutOffsetTime + SpontActivityWindow(2)))))/((SpontActivityWindow(2) - SpontActivityWindow(1))/1000);
        UnDirBoutSpikingDetails.PostBoutWindowSpikeCount(Index) = length(find((ClusterSpikeData >= (BoutOffsetTime + PostSongActivityWindow(1))) & (ClusterSpikeData < (BoutOffsetTime + PostSongActivityWindow(2)))))/((PostSongActivityWindow(2) - PostSongActivityWindow(1))/1000);
        UnDirBoutSpikingDetails.PostMotorBoutWindowSpikeCount(Index) = length(find((ClusterSpikeData >= (BoutOffsetTime + PostMotorSongActivityWindow(1))) & (ClusterSpikeData < (BoutOffsetTime + PostMotorSongActivityWindow(2)))))/((PostMotorSongActivityWindow(2) - PostMotorSongActivityWindow(1))/1000);
        UnDirBoutSpikingDetails.BoutOffsetSpikeTrain{Index} = ClusterSpikeData(find((ClusterSpikeData >= (BoutOffsetTime -200)) & (ClusterSpikeData < (BoutOffsetTime + 2000)))) - BoutOffsetTime;
        UnDirBoutSpikingDetails.BoutSpikeTrain{Index} = ClusterSpikeData(find((ClusterSpikeData >= (BoutOnsetTime - InterBoutInterval)) & (ClusterSpikeData <= (BoutOffsetTime + InterBoutInterval)))) - BoutOffsetTime;
        
        if (~isempty(UnDirBoutSpikingDetails.BoutOffsetSpikeTrain{Index}))
            UnDirBoutSpikingDetails.SpikeRaster{Index} = [UnDirBoutSpikingDetails.BoutOffsetSpikeTrain{Index}(:) ones(size(UnDirBoutSpikingDetails.BoutOffsetSpikeTrain{Index}(:)))];
        else
            UnDirBoutSpikingDetails.SpikeRaster{Index} = [];
        end
        
        UnDirBoutSpikingDetails.BoutLabels{Index} = AllLabels(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        UnDirBoutSpikingDetails.BoutOnsets{Index} = AllOnsets(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        UnDirBoutSpikingDetails.BoutOffsets{Index} = AllOffsets(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        UnDirBoutSpikingDetails.BoutFirstSyllLabel(Index) = UnDirBoutSpikingDetails.BoutLabels{Index}(1);
        UnDirBoutSpikingDetails.BoutPST(Index,:) = conv(histc(ClusterSpikeData, PSTEdges + BoutOffsetTime)/PSTBinSize, GaussWin, 'same');
        UnDirBoutSpikingDetails.SongBoutOrNot(Index) = SURecordingDetails.Bouts(j,7);
        
        Index = Index + 1;
    end
else
    Index = 1;
    for j = UnDirBouts(:)',
        BoutOnsetTime = SURecordingDetails.Bouts(j,5);
        BoutOffsetTime = SURecordingDetails.Bouts(j,6);
 
        UnDirBoutSpikingDetails.BoutIndex(Index) = j;
        
        UnDirBoutSpikingDetails.SpontWindowSpikeCount(Index) = length(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutOffsetTime + SpontActivityWindow(1))) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} < (BoutOffsetTime + SpontActivityWindow(2)))))/((SpontActivityWindow(2) - SpontActivityWindow(1))/1000);
        UnDirBoutSpikingDetails.PostBoutWindowSpikeCount(Index) = length(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutOffsetTime + PostSongActivityWindow(1))) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} < (BoutOffsetTime + PostSongActivityWindow(2)))))/((PostSongActivityWindow(2) - PostSongActivityWindow(1))/1000);
        UnDirBoutSpikingDetails.PostMotorBoutWindowSpikeCount(Index) = length(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutOffsetTime + PostMotorSongActivityWindow(1))) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} < (BoutOffsetTime + PostMotorSongActivityWindow(2)))))/((PostMotorSongActivityWindow(2) - PostMotorSongActivityWindow(1))/1000);
        UnDirBoutSpikingDetails.BoutOffsetSpikeTrain{Index} = SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)}(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutOffsetTime - 200)) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} < (BoutOffsetTime + 2000)))) - BoutOffsetTime;
        UnDirBoutSpikingDetails.BoutSpikeTrain{Index} = SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)}(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutOnsetTime - InterBoutInterval)) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} <= (BoutOffsetTime + InterBoutInterval)))) - BoutOffsetTime;
        
        if (~isempty(UnDirBoutSpikingDetails.BoutOffsetSpikeTrain{Index}))
            UnDirBoutSpikingDetails.SpikeRaster{Index} = [UnDirBoutSpikingDetails.BoutOffsetSpikeTrain{Index}(:) ones(size(UnDirBoutSpikingDetails.BoutOffsetSpikeTrain{Index}(:)))];
        else
            UnDirBoutSpikingDetails.SpikeRaster{Index} = [];
        end
        
        UnDirBoutSpikingDetails.BoutLabels{Index} = SURecordingDetails.NoteInfo{SURecordingDetails.Bouts(j,3)}.labels(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        UnDirBoutSpikingDetails.BoutOnsets{Index} = SURecordingDetails.NoteInfo{SURecordingDetails.Bouts(j,3)}.onsets(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        UnDirBoutSpikingDetails.BoutOffsets{Index} = SURecordingDetails.NoteInfo{SURecordingDetails.Bouts(j,3)}.offsets(SURecordingDetails.Bouts(j,1):SURecordingDetails.Bouts(j,2));
        UnDirBoutSpikingDetails.BoutFirstSyllLabel(Index) = UnDirBoutSpikingDetails.BoutLabels{Index}(1);
        UnDirBoutSpikingDetails.BoutPST(Index,:) = conv(histc(SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)}, PSTEdges + BoutOffsetTime)/PSTBinSize, GaussWin, 'same');
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
        
        if (BoutOffsetTime >= InterBoutInterval)
            DirBoutSpikingDetails.SpontWindowSpikeCount(Index) = length(find((ClusterSpikeData >= (BoutOffsetTime + SpontActivityWindow(1))) & (ClusterSpikeData < (BoutOffsetTime + SpontActivityWindow(2)))))/((SpontActivityWindow(2) - SpontActivityWindow(1))/1000);
            DirBoutSpikingDetails.BoutOffsetSpikeTrain{Index} = ClusterSpikeData(find((ClusterSpikeData >= (BoutOffsetTime - 200)) & (ClusterSpikeData < (BoutOffsetTime + 2000)))) - BoutOffsetTime;
            DirBoutSpikingDetails.BoutPST(Index,:) = conv(histc(ClusterSpikeData, PSTEdges + BoutOnsetTime)/PSTBinSize, GaussWin, 'same');
            if (BoutOffsetTime <= (sum(SURecordingDetails.FileLen) - InterBoutInterval))
                DirBoutSpikingDetails.BoutSpikeTrain{Index} = ClusterSpikeData(find((ClusterSpikeData >= (BoutOnsetTime - InterBoutInterval)) & (ClusterSpikeData <= (BoutOffsetTime + InterBoutInterval)))) - BoutOffsetTime;
            else
                DirBoutSpikingDetails.BoutSpikeTrain{Index} = [];
            end
        else
            DirBoutSpikingDetails.BoutOffsetSpikeTrain{Index} = [];
            DirBoutSpikingDetails.BoutPST(Index,:) = ones(1,length(PSTEdges))*NaN;
            DirBoutSpikingDetails.SpontWindowSpikeCount(Index) = -1; % a negative value to indicate that there was no data
            DirBoutSpikingDetails.BoutSpikeTrain{Index} = [];
        end
        
        if (BoutOffsetTime >= PostSongActivityWindow(1))
            DirBoutSpikingDetails.PostBoutWindowSpikeCount(Index) = length(find((ClusterSpikeData >= (BoutOffsetTime + PostSongActivityWindow(1))) & (ClusterSpikeData < (BoutOffsetTime + PostSongActivityWindow(2)))))/((SpontActivityWindow(2) - SpontActivityWindow(1))/1000);
            UnDirBoutSpikingDetails.PostMotorBoutWindowSpikeCount(Index) = length(find((ClusterSpikeData >= (BoutOffsetTime + PostMotorSongActivityWindow(1))) & (ClusterSpikeData < (BoutOffsetTime + PostMotorSongActivityWindow(2)))))/((PostMotorSongActivityWindow(2) - PostMotorSongActivityWindow(1))/1000);
            DirBoutSpikingDetails.BoutOffsetSpikeTrain{Index} = ClusterSpikeData(find((ClusterSpikeData >= (BoutOffsetTime - 200)) & (ClusterSpikeData < (BoutOffsetTime + 2000)))) - BoutOffsetTime;
            DirBoutSpikingDetails.BoutPST(Index,:) = conv(histc(ClusterSpikeData, PSTEdges + BoutOffsetTime)/PSTBinSize, GaussWin, 'same');
        else
            DirBoutSpikingDetails.PostBoutWindowSpikeCount(Index) = -1;
        end
        
        if (~isempty(DirBoutSpikingDetails.BoutOffsetSpikeTrain{Index}))
            DirBoutSpikingDetails.SpikeRaster{Index} = [DirBoutSpikingDetails.BoutOffsetSpikeTrain{Index}(:) ones(size(DirBoutSpikingDetails.BoutOffsetSpikeTrain{Index}(:)))];
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
        
        if (BoutOffsetTime >= InterBoutInterval)
            DirBoutSpikingDetails.SpontWindowSpikeCount(Index) = length(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutOffsetTime + SpontActivityWindow(1))) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} < (BoutOffsetTime + SpontActivityWindow(2)))))/((SpontActivityWindow(2) - SpontActivityWindow(1))/1000);
            if ((BoutOffsetTime + 200) <= SURecordingDetails.FileLen(SURecordingDetails.Bouts(j,3)))
                DirBoutSpikingDetails.BoutOffsetSpikeTrain{Index} = SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)}(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutOffsetTime - 200)) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} < (BoutOffsetTime + 2000)))) - BoutOffsetTime;
                DirBoutSpikingDetails.BoutPST(Index,:) = conv(histc(SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)}, PSTEdges + BoutOffsetTime)/PSTBinSize, GaussWin, 'same');
            else
                DirBoutSpikingDetails.BoutOffsetSpikeTrain{Index} = [];
                DirBoutSpikingDetails.BoutPST(Index,:) = ones(1, length(PSTEdges))*NaN;
            end
            if ((BoutOffsetTime + InterBoutInterval) <= SURecordingDetails.FileLen(SURecordingDetails.Bouts(j,3)))
                DirBoutSpikingDetails.BoutSpikeTrain{Index} = SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)}(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutOnsetTime - InterBoutInterval)) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} <= (BoutOffsetTime + InterBoutInterval)))) - BoutOffsetTime;
            else
                DirBoutSpikingDetails.BoutSpikeTrain{Index} = [];
            end
        else
            DirBoutSpikingDetails.BoutOffsetSpikeTrain{Index} = [];
            DirBoutSpikingDetails.BoutPST(Index,:) = ones(1,length(PSTEdges))*NaN;
            DirBoutSpikingDetails.SpontWindowSpikeCount(Index) = -1; % a negative value to indicate that there was no data
            DirBoutSpikingDetails.BoutSpikeTrain{Index} = [];
        end
        
        if (BoutOffsetTime >= PostSongActivityWindow(1))
            DirBoutSpikingDetails.PostBoutWindowSpikeCount(Index) = length(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutOffsetTime + PostSongActivityWindow(1))) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} < (BoutOffsetTime + PostSongActivityWindow(2)))))/((SpontActivityWindow(2) - SpontActivityWindow(1))/1000);
            DirBoutSpikingDetails.PostMotorBoutWindowSpikeCount(Index) = length(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutOffsetTime + PostMotorSongActivityWindow(1))) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} < (BoutOffsetTime + PostMotorSongActivityWindow(2)))))/((PostMotorSongActivityWindow(2) - PostMotorSongActivityWindow(1))/1000);
            DirBoutSpikingDetails.BoutOffsetSpikeTrain{Index} = SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)}(find((SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} >= (BoutOffsetTime - 200)) & (SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)} < (BoutOffsetTime + 2000)))) - BoutOffsetTime;
            DirBoutSpikingDetails.BoutPST(Index,:) = conv(histc(SURecordingDetails.ClusterSpikeTimes{SURecordingDetails.Bouts(j,3)}, PSTEdges + BoutOffsetTime)/PSTBinSize, GaussWin, 'same');
            % In this case, bout onset spike train will have only the times
            % for the pre-song window. If spont. window spike count is -1,
            % then I know that the raster should be limited to the pre-song
            % window alone
        else
            DirBoutSpikingDetails.PostBoutWindowSpikeCount(Index) = -1;
        end
        
        if (~isempty(DirBoutSpikingDetails.BoutOffsetSpikeTrain{Index}))
            DirBoutSpikingDetails.SpikeRaster{Index} = [DirBoutSpikingDetails.BoutOffsetSpikeTrain{Index}(:) ones(size(DirBoutSpikingDetails.BoutOffsetSpikeTrain{Index}(:)))];
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