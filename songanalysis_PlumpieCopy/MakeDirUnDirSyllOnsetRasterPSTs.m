function [SyllOnsetRasterPST] = MakeDirUnDirSyllOnsetRasterPSTs(BoutDetails, SURecordingDetails)

% First find all the motif syllabless in each bout and see if there is enough pretime and post time
PrePostMotifTime = 50; % in ms

Fs = 10000;
GaussianLen = 2;
Width = 0.005; % in sec

XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * Fs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs) * (Width * Fs)));

% For the smoothed PST
Time = -PrePostMotifTime/1000:1/Fs:PrePostMotifTime/1000;

for i = 1:length(SURecordingDetails.MotifLabels),
    SyllOnsetRasterPST.Time{i} = Time;
    ValidSongBouts = find(~isnan(BoutDetails.DirUnDir));
    SyllOnsetRasterPST.SyllLabel(i) = SURecordingDetails.MotifLabels(i);
    for j = ValidSongBouts(:)',
        Motifs = find(BoutDetails.BoutLabels{j} == SURecordingDetails.MotifLabels(i));
        for k = 1:length(Motifs),
            if ((BoutDetails.BoutOnsets{j}(Motifs(k)) >= PrePostMotifTime) && ((BoutDetails.BoutOnsets{j}(Motifs(k)) + PrePostMotifTime) <= BoutDetails.FileEndTime(j)))
                SyllOnsetRasterPST.MotifStartIndices{i}{j}(k) = Motifs(k);
                SyllOnsetRasterPST.Onsets{i}{j}(k,:) = BoutDetails.BoutOnsets{j}(Motifs(k)) - BoutDetails.BoutOnsets{j}(Motifs(k));
                SyllOnsetRasterPST.Offsets{i}{j}(k,:) = BoutDetails.BoutOffsets{j}(Motifs(k)) - BoutDetails.BoutOnsets{j}(Motifs(k));
                SyllOnsetRasterPST.Raster{i}{j}{k} = BoutDetails.SpikeTimes{j}(find((BoutDetails.SpikeTimes{j} >= (BoutDetails.BoutOnsets{j}(Motifs(k)) - PrePostMotifTime)) & (BoutDetails.SpikeTimes{j} <= (BoutDetails.BoutOnsets{j}(Motifs(k)) + PrePostMotifTime))));
                SyllOnsetRasterPST.Raster{i}{j}{k} = SyllOnsetRasterPST.Raster{i}{j}{k} - BoutDetails.BoutOnsets{j}(Motifs(k));
                % Now for the smoothed PST
                FR = zeros(1,length(Time));
                for Spike = 1:length(SyllOnsetRasterPST.Raster{i}{j}{k}),
                    if (ceil((PrePostMotifTime + SyllOnsetRasterPST.Raster{i}{j}{k}(Spike)) * Fs/1000) > 0)
                        FR(1,ceil((PrePostMotifTime + SyllOnsetRasterPST.Raster{i}{j}{k}(Spike)) * Fs/1000)) = 1;
                    else
                        FR(1,1) = 1;
                    end
                end
                SyllOnsetRasterPST.SmoothedPST{i}{j}{k} = conv(FR, GaussWin, 'same');
            end
        end
    end
    if (isempty(ValidSongBouts))
        SyllOnsetRasterPST.MotifStartIndices{i} = [];
        SyllOnsetRasterPST.Onsets{i} = [];
        SyllOnsetRasterPST.Offsets{i} = [];
        SyllOnsetRasterPST.Raster{i} = [];
        SyllOnsetRasterPST.Time{i} = [];
        SyllOnsetRasterPST.SmoothedPST{i} = [];
    end
end