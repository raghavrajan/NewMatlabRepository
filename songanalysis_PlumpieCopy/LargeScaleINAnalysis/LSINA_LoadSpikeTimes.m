function [SpikeInfo] = LSINA_LoadSpikeTimes(BirdParameters)

% Now for each of the filelists, load up the spiketimes
for j = 1:length(BirdParameters.SongFileNames),
    DotIndex = find(BirdParameters.SongFileList == '.');
    SpikeFileDir = BirdParameters.SongFileList(1:DotIndex(end)-1);
    SpikeInfo{j} = load(fullfile(BirdParameters.DataDirectory, SpikeFileDir, ['Chan', num2str(BirdParameters.SpikeChanNo), '_', BirdParameters.SongFileNames{j}, '.spk']));
    % Check to see if the times are in ms or sec
    % Right now, I am using a hack that if all the SpikeTimes are < 1000,
    % then it must be in sec and so it should then be converted into ms by
    % multiplying by 1000
    if (~isempty(SpikeInfo{j}))
        if (isempty(find(SpikeInfo{j}(:,2) > 1000, 1, 'first')))
            SpikeInfo{j}(:,2) = SpikeInfo{j}(:,2)*1000;
        end
    end
end