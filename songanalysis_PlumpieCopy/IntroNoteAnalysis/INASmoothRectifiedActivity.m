function [SpikeData, ASpikeData] = INASmoothRectifiedActivity(DataDir, RecFileDir, NoteFiles, FileType, GaussWin, ContinuousOrNot, SpikeChanNo)

PresentDir = pwd;

if (strfind(ContinuousOrNot, 'continuous'))
    CSpikeData = [];
else
    SpikeData = [];
end

for i = 1:length(NoteFiles),
    Index = strfind(NoteFiles(i).name, '.not.mat');
    SongFile = NoteFiles(i).name(1:Index-1);
    disp(SongFile);
    
    if (strfind(FileType, 'okrank'))
        [Data, Fs] = SSAReadOKrankData(DataDir, RecFileDir, SongFile, SpikeChanNo);
    else
        if (strfind(FileType, 'obs'));
            [Data, Fs] = SSASoundIn(DataDir, RecFileDir, SongFile, ['obs', num2str(SpikeChanNo), 'r']);
            Data = Data * 5/32768;
        end
    end
    SpikeData{i} = conv(Data.^2, GaussWin, 'same');
end
ASpikeData = SpikeData;

figure
subplot(2,1,1);
plot(Data);
subplot(2,1,2);
plot(SpikeData{end});
if (strfind(ContinuousOrNot, 'continuous'))
    SpikeData = CSpikeData;
end

cd(PresentDir);