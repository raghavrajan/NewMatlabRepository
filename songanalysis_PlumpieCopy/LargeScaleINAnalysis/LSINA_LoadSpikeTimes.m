function [SpikeInfo] = LSINA_LoadSpikeTimes(BirdParameters)

% Now for each of the filelists, load up the spiketimes
for j = 1:length(BirdParameters.SongFileNames),
    DotIndex = find(BirdParameters.SongFileList == '.');
    SpikeFileDir = BirdParameters.SongFileList(1:DotIndex(end)-1);
    SpikeInfo{j} = load(fullfile(BirdParameters.DataDirectory, SpikeFileDir, ['Chan', num2str(BirdParameters.SpikeChanNo), '_', BirdParameters.SongFileNames{j}, '.spk']));
end