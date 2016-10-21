function [UnDirSpikeTrain, UnDirIFR, UnDirPST, UnDirRaster, DirSpikeTrain, DirIFR, DirPST, DirRaster] = LSINA_MakeRasterPST_NonSongBouts(SURecordingDetails, NeuronIndex, InterboutInterval, PSTEdges, IFRTime, GaussWin)

SongBouts_WithPreData = find((SURecordingDetails(NeuronIndex).Bouts(:,7) == 0) & (SURecordingDetails(NeuronIndex).Bouts(:,8) >= 0));
% For each of these bouts check if it is directed or not and then
% accordingly pull out the spikes
DirBoutNo = 0;
UnDirBoutNo = 0;
DirRaster = [];
UnDirRaster = [];
DirSpikeTrain = [];
UnDirSpikeTrain = [];
DirPST = [];
UnDirPST = [];
DirIFR = [];
UnDirIFR = [];


for i = 1:length(SongBouts_WithPreData),

    BoutOnsetFile = SURecordingDetails(NeuronIndex).Bouts(SongBouts_WithPreData(i),3);
    % First check if this is directed or not
    DirUndir(i) = 0; % 0 for undir and 1 for dir
    for j = 1:size(SURecordingDetails(NeuronIndex).DirPresentations,1),
        if ((BoutOnsetFile >= SURecordingDetails(NeuronIndex).DirPresentations(j,1)) & (BoutOnsetFile <= SURecordingDetails(NeuronIndex).DirPresentations(j,4)))
            DirUndir(i) = 1;
            break;
        end
    end

    % Bout Onset refers to the first syllable of the bout. Now, I have
    % to check if 2000ms before this onset is still in the same file or
    % is now in the previous file
    BoutOnset = SURecordingDetails(NeuronIndex).Bouts(SongBouts_WithPreData(i),5);


    if (BoutOnset >= InterboutInterval)
        % Now to check if there is atleast 200ms after the start of the
        % bout in the same file
        if ((BoutOnset + 200) <= SURecordingDetails(NeuronIndex).FileLen(BoutOnsetFile))
            TrialSpikeIndices = find((SURecordingDetails(NeuronIndex).ClusterSpikeTimes{BoutOnsetFile} >= (BoutOnset - InterboutInterval)) & (SURecordingDetails(NeuronIndex).ClusterSpikeTimes{BoutOnsetFile} <= (BoutOnset + 200)));
            if (~isempty(TrialSpikeIndices))
                TrialSpikeTrain = SURecordingDetails(NeuronIndex).ClusterSpikeTimes{BoutOnsetFile}(TrialSpikeIndices) - BoutOnset;
            else
                TrialSpikeTrain = [];
            end
        else
            if (BoutOnsetFile == length(SURecordingDetails(NeuronIndex).SongFileNames))
                continue;
            end
            % This means there is not enough post data in the same
            % file. So then get all the spikes from this file from
            % (BoutOnset - 2000) to end. Then get the extra data from
            % the next file
            TrialSpikeIndices = find(SURecordingDetails(NeuronIndex).ClusterSpikeTimes{BoutOnsetFile} >= (BoutOnset - 2000));
            if (~isempty(TrialSpikeIndices))
                TrialSpikeTrain = SURecordingDetails(NeuronIndex).ClusterSpikeTimes{BoutOnsetFile}(TrialSpikeIndices) - BoutOnset;
            else
                TrialSpikeTrain = [];
            end
            % Now get the extra spikes
            ExtraTime = (BoutOnset + 200) - SURecordingDetails(NeuronIndex).FileLen(BoutOnsetFile);
            TrialSpikeIndices = find(SURecordingDetails(NeuronIndex).ClusterSpikeTimes{BoutOnsetFile + 1} <= ExtraTime);
            if (~isempty(TrialSpikeIndices))
                ExtraSpikeTrain = SURecordingDetails(NeuronIndex).ClusterSpikeTimes{BoutOnsetFile + 1}(TrialSpikeIndices) + SURecordingDetails(NeuronIndex).FileLen(BoutOnsetFile) - BoutOnset;
            else
                ExtraSpikeTrain = [];
            end
            TrialSpikeTrain = [TrialSpikeTrain(:); ExtraSpikeTrain(:)];
            clear ExtraSpikeTrain;
        end
    else
        % This means I need to get some pre data from the previous file
        PreExtraTime = InterboutInterval - BoutOnset;
        TrialSpikeIndices = find(SURecordingDetails(NeuronIndex).ClusterSpikeTimes{BoutOnsetFile - 1} >= (SURecordingDetails(NeuronIndex).FileLen(BoutOnsetFile - 1) - PreExtraTime));
        if (~isempty(TrialSpikeIndices))
            PreExtraSpikeTrain = -(SURecordingDetails(NeuronIndex).ClusterSpikeTimes{BoutOnsetFile - 1}(TrialSpikeIndices)) + SURecordingDetails(NeuronIndex).FileLen(BoutOnsetFile - 1) - PreExtraTime;
        else
            PreExtraSpikeTrain = [];
        end

        % Now see if the offset is fine or not.
        if ((BoutOnset + 200) <= SURecordingDetails(NeuronIndex).FileLen(BoutOnsetFile))
            TrialSpikeIndices = find((SURecordingDetails(NeuronIndex).ClusterSpikeTimes{BoutOnsetFile} <= (BoutOnset + 200)));
            if (~isempty(TrialSpikeIndices))
                TrialSpikeTrain = SURecordingDetails(NeuronIndex).ClusterSpikeTimes{BoutOnsetFile}(TrialSpikeIndices) - BoutOnset;
            else
                TrialSpikeTrain = [];
            end
        else
            % This means there is not enough post data in the same
            % file. So then get all the spikes from this file from
            % (BoutOnset - 2000) to end. Then get the extra data from
            % the next file
            TrialSpikeIndices = 1:1:length(SURecordingDetails(NeuronIndex).ClusterSpikeTimes{BoutOnsetFile});
            if (~isempty(TrialSpikeIndices))
                TrialSpikeTrain = SURecordingDetails(NeuronIndex).ClusterSpikeTimes{BoutOnsetFile}(TrialSpikeIndices) - BoutOnset;
            else
                TrialSpikeTrain = [];
            end
            % Now get the extra spikes
            ExtraTime = (BoutOnset + 200) - SURecordingDetails(NeuronIndex).FileLen(BoutOnsetFile);
            TrialSpikeIndices = find(SURecordingDetails(NeuronIndex).ClusterSpikeTimes{BoutOnsetFile + 1} <= ExtraTime);
            if (~isemty(TrialSpikeIndices))
                ExtraSpikeTrain = SURecordingDetails(NeuronIndex).ClusterSpikeTimes{BoutOnsetFile + 1}(TrialSpikeIndices) + SURecordingDetails(NeuronIndex).FileLen(BoutOnsetFile) - BoutOnset;
            else
                ExtraSpikeTrain = [];
            end
            TrialSpikeTrain = [TrialSpikeTrain(:); ExtraSpikeTrain(:)];
        end
        TrialSpikeTrain = [PreExtraSpikeTrain(:); TrialSpikeTrain(:)];
        clear PreExtraSpikeTrain;
    end

    TrialSpikeTrain = TrialSpikeTrain/1000; 

    if (DirUndir(i) == 0)
        UnDirBoutNo = UnDirBoutNo + 1;
        UnDirSpikeTrain{UnDirBoutNo} = TrialSpikeTrain;
        UnDirIFR(UnDirBoutNo,:) = conv(CalculateIFR(TrialSpikeTrain + (InterboutInterval/1000), IFRTime), GaussWin, 'same');
        if (isempty(TrialSpikeTrain))
            UnDirPST(UnDirBoutNo,:) = zeros(size(PSTEdges(:)'));
        else
            UnDirPST(UnDirBoutNo,:) = histc(TrialSpikeTrain, PSTEdges);
        end
        UnDirRaster = [UnDirRaster; [TrialSpikeTrain(:) ones(size(TrialSpikeTrain(:)))*UnDirBoutNo]];
    else
        DirBoutNo = DirBoutNo + 1;
        DirSpikeTrain{DirBoutNo} = TrialSpikeTrain;
        DirRaster = [DirRaster; [TrialSpikeTrain(:) ones(size(TrialSpikeTrain(:)))*DirBoutNo]];
        if (isempty(TrialSpikeTrain))
            DirPST(DirBoutNo,:) = zeros(size(PSTEdges(:)'));
        else
            DirPST(DirBoutNo,:) = histc(TrialSpikeTrain, PSTEdges);
        end
        DirIFR(DirBoutNo,:) = conv(CalculateIFR(TrialSpikeTrain  + (InterboutInterval/1000), IFRTime), GaussWin, 'same');
    end
end