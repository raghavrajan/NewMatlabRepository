function [WarpedRasterPSTs] = MakeWarpedDirUnDirRasterPSTs(CombinedData)

PSTBinSize = 3;

Fs = 10000;
GaussianLen = 2;
Width = 0.005; % in sec

XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * Fs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs) * (Width * Fs)));


for i = 1:length(CombinedData.Offsets),
    WarpedRasterPSTs.WarpedRaster{i} = [];
    % First find median motif
    [SortedDurs, SortedIndices] = sort(CombinedData.Offsets{i}(:,end));
    MedianMotifIndex = SortedIndices(round(length(SortedIndices)/2));
    
    MedianMotifOnsetsOffsets = [CombinedData.Onsets{i}(MedianMotifIndex,:); CombinedData.Offsets{i}(MedianMotifIndex,:)];

    % For the smoothed PST
    Time = -0.1:1/Fs:MedianMotifOnsetsOffsets(2,end)/1000;

    PSTEdges = -100:PSTBinSize:CombinedData.Offsets{i}(MedianMotifIndex, end);
    WarpedRasterPSTs.PSTEdges{i} = PSTEdges;
    WarpedRasterPSTs.MedianMotifOnsetsOffsets{i} = MedianMotifOnsetsOffsets;
    WarpedRasterPSTs.Time{i} = Time;
    
    for j = 1:size(CombinedData.Onsets{i},1),
        WarpedRasterPSTs.WarpedRaster{i}{j} = [];
        % Now for warping, first find the spikes before 0 and don't warp
        % these
        Spikes = CombinedData.Raster{i}{j}(find(CombinedData.Raster{i}{j} < 0));
        WarpedRasterPSTs.WarpedRaster{i}{j} = [WarpedRasterPSTs.WarpedRaster{i}{j}; Spikes(:)];
        
        % Now find the spikes in each syllable and the interval after and
        % warp them to the median interval
        for Syll = 1:size(CombinedData.Onsets{i},2),
            Spikes = CombinedData.Raster{i}{j}(find((CombinedData.Raster{i}{j} >= CombinedData.Onsets{i}(j,Syll)) & (CombinedData.Raster{i}{j} < CombinedData.Offsets{i}(j,Syll))));
            if (~isempty(Spikes))
                % First make it on a scale of 0 to 1 based on the duration of
                % the syllable it is from
                Spikes = ((Spikes - CombinedData.Onsets{i}(j,Syll))/(CombinedData.Offsets{i}(j,Syll) - CombinedData.Onsets{i}(j,Syll))); 
                % Now put it onto the timescale of the median motif syllable
                % duration given above
                Spikes = MedianMotifOnsetsOffsets(1,Syll) + (Spikes * diff(MedianMotifOnsetsOffsets(:,Syll)));
                WarpedRasterPSTs.WarpedRaster{i}{j} = [WarpedRasterPSTs.WarpedRaster{i}{j}; Spikes(:)];
            end            
            % Now do the same for the interval between the syllables for
            % all except the last
            if (Syll < size(CombinedData.Onsets{i},2))
                Spikes = CombinedData.Raster{i}{j}(find((CombinedData.Raster{i}{j} >= CombinedData.Offsets{i}(j,Syll)) & (CombinedData.Raster{i}{j} < CombinedData.Onsets{i}(j,Syll+1))));
                if (~isempty(Spikes))
                    Spikes = ((Spikes - CombinedData.Offsets{i}(j,Syll))/(CombinedData.Onsets{i}(j,Syll+1) - CombinedData.Offsets{i}(j,Syll))); 
                    Spikes = MedianMotifOnsetsOffsets(2,Syll) + (Spikes * (MedianMotifOnsetsOffsets(1,Syll+1) - MedianMotifOnsetsOffsets(2,Syll)));
                    WarpedRasterPSTs.WarpedRaster{i}{j} = [WarpedRasterPSTs.WarpedRaster{i}{j}; Spikes(:)];
                end
            end
        end
        if (~isempty(WarpedRasterPSTs.WarpedRaster{i}{j}))
            WarpedRasterPSTs.PST{i}(j,:) = histc(WarpedRasterPSTs.WarpedRaster{i}{j}, PSTEdges)/(PSTBinSize/1000);
        else
            WarpedRasterPSTs.PST{i}(j,:) = zeros(1, length(PSTEdges));
        end
        % Now for the smoothed PST
        FR = zeros(1,length(Time));
        for Spike = 1:length(WarpedRasterPSTs.WarpedRaster{i}{j}),
            if (ceil((100 + WarpedRasterPSTs.WarpedRaster{i}{j}(Spike)) * Fs/1000) > 0)
                FR(1,ceil((100 + WarpedRasterPSTs.WarpedRaster{i}{j}(Spike)) * Fs/1000)) = 1;
            else
                FR(1,1) = 1;
            end
        end
        WarpedRasterPSTs.SmoothedPST{i}(j,:) = conv(FR, GaussWin, 'same');
    end
end