function [AverageCorr] = MA_DoSpectralCorrelations(Bouts, SongFileNames, RawDataDir, FileType, WindowSize, StepSize, BandWidth, BoutSize)

Bouts = cell2mat(Bouts);
BoutLens = Bouts(:,2) - Bouts(:,1);

LongBouts = find(BoutLens >= (BoutSize*1000));

if (length(LongBouts) <= 10)
    ActualBouts = LongBouts;
else
    RandomIndices = randperm(length(LongBouts));
    ActualBouts = LongBouts(RandomIndices(1:10));
end

disp(['    Calculating multi-taper power spectral densities for ', num2str(length(ActualBouts)), ' random bouts ...']);

for i = 1:length(ActualBouts),
    [Song, Fs] = MA_ReadSongFile(RawDataDir, SongFileNames{Bouts(ActualBouts(i),4)}, FileType);
    BoutRawData{i} = Song(ceil(Bouts(ActualBouts(i),1)*Fs/1000):floor(Bouts(ActualBouts(i),2)*Fs/1000));
    [PowSpect{i}, Freq{i}] = CalculateMultiTaperSpectrogram(BoutRawData{i}, Fs, WindowSize, StepSize, BandWidth);
    PowSpect{i} = PowSpect{i}(find((Freq{i} >= 860) & (Freq{i} <= 8600)),:);
end

disp(['    Calculating spectral correlations ...']);
Index = 0;
for i = 1:10,
    for j = i+1:10,
        SpectralCorrs{i}{j} = corr(PowSpect{i}, PowSpect{j});
        Index = Index + 1;
        for k = -1500:1:1500,
            AverageCorr(Index, k+1501) = mean(diag(SpectralCorrs{i}{j}, k));
        end
    end
end
