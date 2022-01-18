function [MeanMotifStartActivity] = PlotGaussianSmoothedNeuralActivity(DataDir, FileList, FileType, SpikeChanNo, Motif, INs)

Fid = fopen(FileList, 'r');
Files = textscan(Fid, '%s', 'DeLimiter', '\n');
Files = Files{1};
fclose(Fid);

[RawData, Fs] = GetData(DataDir, Files{1}, FileType, 0);

Width = 0.002; % in seconds - Gaussian width or std
GaussianLen = 3; % length of gaussian

XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * Fs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs) * (Width * Fs)));

for i = 1:length(Files),
    NoteInfo{i} = load(fullfile(DataDir, 'ASSLNoteFiles', [Files{i}, '.not.mat']));
    for j = 1:length(SpikeChanNo),
        [SpikeData, Fs] = GetData(DataDir, Files{i}, FileType, SpikeChanNo(j));
        RawSpikeData{i}{j} = SpikeData;
        SmoothedSpikeData{i}{j} = conv(abs(SpikeData), GaussWin, 'same');
    end
    [SongData, Fs] = GetData(DataDir, Files{i}, FileType, 1);
    SmoothedSongData{i}{j} = conv(abs(SongData), GaussWin, 'same');
end

% Now find all the files with motifs that begin after 2s and don't have
% calls before them and get the times of motif start. Then get out the
% smoothed neural data starting from 2s before motif start to 1s after
SongStartActivity = [];
for i = 1:length(SpikeChanNo),
    MotifStartActivity{i} = [];
    RawSpikeMotifStartActivity{i} = [];
end
PreTime = 2; % in seconds
PostTime = 1; % in seconds

for i = 1:length(NoteInfo),
    Motifs = strfind(NoteInfo{i}.labels, Motif);
    if (~isempty(Motifs))
        % Now check that there are no syllables before the first syllable
        % of the motif and that there is > 2s of data before
        if (isempty(INs))
            if ((Motifs(1) == 1) && (NoteInfo{i}.onsets(Motifs(1))/1000 > PreTime))
                MotifOnsetTime = NoteInfo{i}.onsets(Motifs(1))/1000;
                for j = 1:length(SpikeChanNo),
                    MotifStartActivity{j}(end+1,:) = SmoothedSpikeData{i}{j}(round((MotifOnsetTime - PreTime) * Fs):round((MotifOnsetTime + PostTime) * Fs));
                    RawSpikeMotifStartActivity{j}(end+1,:) = RawSpikeData{i}{j}(round((MotifOnsetTime - PreTime) * Fs):round((MotifOnsetTime + PostTime) * Fs));
                end
                SongStartActivity(end+1,:) = SmoothedSongData{i}{j}(round((MotifOnsetTime - PreTime) * Fs):round((MotifOnsetTime + PostTime) * Fs));
            end
        else
            if ((length(find(NoteInfo{i}.labels(1:Motifs(1)) == INs)) == (Motifs(1) - 1)) && (NoteInfo{i}.onsets(1)/1000 > PreTime))
                MotifOnsetTime = NoteInfo{i}.onsets(Motifs(1))/1000;
                for j = 1:length(SpikeChanNo),
                    MotifStartActivity{j}(end+1,:) = SmoothedSpikeData{i}{j}(round((MotifOnsetTime - PreTime) * Fs):round((MotifOnsetTime + PostTime) * Fs));
                    RawSpikeMotifStartActivity{j}(end+1,:) = RawSpikeData{i}{j}(round((MotifOnsetTime - PreTime) * Fs):round((MotifOnsetTime + PostTime) * Fs));
                end
                SongStartActivity(end+1,:) = SmoothedSongData{i}{j}(round((MotifOnsetTime - PreTime) * Fs):round((MotifOnsetTime + PostTime) * Fs));
            end
        end
    end
end

PSTHTime = linspace(-PreTime, PostTime, size(MotifStartActivity{1}, 2));
Colours = distinguishable_colors(length(SpikeChanNo)+1);
figure;
hold on;
for i = 1:length(SpikeChanNo),
    if (i == 1)
        [Ax, H1, H2] = plotyy(PSTHTime, mean(SongStartActivity), PSTHTime, mean(MotifStartActivity{i}));
        set(H1, 'Color', 'w');
        set(H2, 'Color', 'w');
        axes(Ax(2));
        hold on;
        patch([PSTHTime fliplr(PSTHTime)], [(mean(MotifStartActivity{i}) + (std(MotifStartActivity{i})/sqrt(size(MotifStartActivity{i}, 1)))) fliplr(mean(MotifStartActivity{i}) - (std(MotifStartActivity{i})/sqrt(size(MotifStartActivity{i}, 1))))], 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.1, 'FaceColor', Colours(i+1,:));
        plot(PSTHTime, mean(MotifStartActivity{i}), 'k', 'Color', Colours(i+1,:));
        axes(Ax(1));
        hold on;
        plot(PSTHTime, mean(SongStartActivity), 'k', 'Color', Colours(1,:));
    else
        axes(Ax(2));
        plot(PSTHTime, mean(MotifStartActivity{i}), 'k', 'Color', Colours(i+1,:));
        patch([PSTHTime fliplr(PSTHTime)], [(mean(MotifStartActivity{i}) + (std(MotifStartActivity{i})/sqrt(size(MotifStartActivity{i}, 1)))) fliplr(mean(MotifStartActivity{i}) - (std(MotifStartActivity{i})/sqrt(size(MotifStartActivity{i}, 1))))], 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.1, 'FaceColor', Colours(i+1,:));
        plot(PSTHTime, mean(MotifStartActivity{i}), 'k', 'Color', Colours(i+1,:));
    end
end
axes(Ax(1));
axis tight;
Temp = axis;
Temp = [-PreTime PostTime 0 Temp(4)*1.05];
axis(Temp);

axes(Ax(2));
axis tight;
Temp = axis;
Temp = [-PreTime PostTime 0 Temp(4)*1.05];
axis(Temp);
patch([0 0 PostTime PostTime], [0 Temp(4) Temp(4) 0], 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.1);

figure;
p = panel();
p.pack(length(SpikeChanNo)+1, 1);

p(1,1).select();
imagesc(PSTHTime, 1:1:size(SongStartActivity, 1), SongStartActivity);
axis([-PreTime PostTime 0.5 (size(SongStartActivity,1)+0.5)]);

for i = 1:length(SpikeChanNo),
    p(i+1, 1).select();
    imagesc(PSTHTime, 1:1:size(MotifStartActivity{i}, 1), MotifStartActivity{i});
    axis([-PreTime PostTime 0.5 (size(MotifStartActivity{i},1)+0.5)]);
end

MeanMotifStartActivity = mean(MotifStartActivity{1});
disp('Finished');