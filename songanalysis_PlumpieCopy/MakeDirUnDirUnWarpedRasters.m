function [RasterPST] = MakeDirUnDirUnWarpedRasters(BoutDetails, SURecordingDetails)

% First find all the motifs in each bout and see if there is enough pretime and post time
PrePostMotifTime = 100; % in ms

for i = 1:length(SURecordingDetails.CommonMotifs),
    ValidSongBouts = find(~isnan(BoutDetails.DirUnDir));
    for j = ValidSongBouts(:)',
        Motifs = strfind(BoutDetails.BoutLabels{j}(:)', SURecordingDetails.CommonMotifs{i});
        for k = 1:length(Motifs),
            if ((BoutDetails.BoutOnsets{j}(Motifs(k)) >= PrePostMotifTime) && ((BoutDetails.BoutOffsets{j}(Motifs(k) + length(SURecordingDetails.CommonMotifs{i}) - 1) + PrePostMotifTime) <= BoutDetails.FileEndTime(j)))
                RasterPST.MotifStartIndices{i}{j}(k) = Motifs(k);
                RasterPST.Onsets{i}{j}(k,:) = BoutDetails.BoutOnsets{j}(Motifs(k):Motifs(k)+length(SURecordingDetails.CommonMotifs{i})-1) - BoutDetails.BoutOnsets{j}(Motifs(k));
                RasterPST.Offsets{i}{j}(k,:) = BoutDetails.BoutOffsets{j}(Motifs(k):Motifs(k)+length(SURecordingDetails.CommonMotifs{i})-1) - BoutDetails.BoutOnsets{j}(Motifs(k));
                RasterPST.Raster{i}{j}{k} = BoutDetails.SpikeTimes{j}(find((BoutDetails.SpikeTimes{j} >= (BoutDetails.BoutOnsets{j}(Motifs(k)) - PrePostMotifTime)) & (BoutDetails.SpikeTimes{j} <= (BoutDetails.BoutOffsets{j}(Motifs(k) + length(SURecordingDetails.CommonMotifs{i}) - 1)))));
                RasterPST.Raster{i}{j}{k} = RasterPST.Raster{i}{j}{k} - BoutDetails.BoutOnsets{j}(Motifs(k));
            end
        end
    end
    if (isempty(ValidSongBouts))
        RasterPST.MotifStartIndices{i} = [];
        RasterPST.Onsets{i} = [];
        RasterPST.Offsets{i} = [];
        RasterPST.Raster{i} = [];
    end
end