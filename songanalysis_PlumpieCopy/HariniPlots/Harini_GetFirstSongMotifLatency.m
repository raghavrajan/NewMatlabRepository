function [] = Harini_GetCumulativeBoutNos(CSVTextFile)

ContinuousFileTime = 30; % in seconds

% First pre-processing step
[BirdParameters, Flag] = Harini_ProcessSongData(CSVTextFile, 2000);

% Now I need to get the first valid song bout and calculate the latency to
% the first song

% I also need the total number of song bouts. In this case, there are song
% bouts that sometimes spread over two files. In this case I need to check
% the time of the file and see whether two files are continuous or not.

% For this, I need to first get the times for each of the files
for i = 1:length(BirdParameters),
    for j = 1:length(BirdParameters(i).SongFileNames),
        ExtensionIndex = strfind(BirdParameters(i).SongFileNames{j}, '.wav');
        FileTimeString = BirdParameters(i).SongFileNames{j}((ExtensionIndex - 6):(ExtensionIndex - 1)); % string from filename in the format hhmmss (h for hour, m for minutes, s for seconds)
        BirdParameters(i).FileTime(j) = (str2double(FileTimeString(1:2))) + (str2double(FileTimeString(3:4))/60) + (str2double(FileTimeString(5:6))/3600); % in hours
    end
end

for i = 1:length(BirdParameters),
    FemaleIntroductionTime = (str2double(BirdParameters(i).Femaleintroductiontime(1:2))) + (str2double(BirdParameters(i).Femaleintroductiontime(3:4))/60) + (str2double(BirdParameters(i).Femaleintroductiontime(5:6))/3600); % in hours
    FemaleExitTime = (str2double(BirdParameters(i).Femaleexittime(1:2))) + (str2double(BirdParameters(i).Femaleexittime(3:4))/60) + (str2double(BirdParameters(i).Femaleexittime(5:6))/3600); % in hours
    SongBouts = find(BirdParameters(i).Bouts(:,7) == 1);
    if (~isempty(SongBouts))
        FirstSongBoutTime = BirdParameters(i).FileTime(BirdParameters(i).Bouts(SongBouts(1),3)) + (BirdParameters(i).Bouts(SongBouts(1),5)/(3600*1000)); % in hours
        for j = BirdParameters(i).Bouts(SongBouts(1),1):BirdParameters(i).Bouts(SongBouts(1),2),
            if (~isempty(find(BirdParameters(i).MotifLabels == char(BirdParameters(i).SyllableData(j,1)))))
                FirstMotifTime = BirdParameters(i).FileTime(BirdParameters(i).Bouts(SongBouts(1),3)) + (BirdParameters(i).SyllableData(j,4)/(3600*1000)); % in hours
                break;
            else
                FirstMotifTime = NaN;
            end
        end
        LatencyToFirstMotif(i) = (FirstMotifTime - FemaleIntroductionTime) * 60; % in minutes
        LatencyToFirstSong(i) = (FirstSongBoutTime - FemaleIntroductionTime) * 60; % in minutes
    else
        LatencyToFirstSong(i) = (FemaleExitTime - FemaleIntroductionTime) * 60; % in minutes
        LatencyToFirstMotif(i) = (FemaleExitTime - FemaleIntroductionTime) * 60; % in minutes
    end
    if (LatencyToFirstSong(i) < 0)
        disp([BirdParameters(i).SongFileNames{1}, ': female intro time = ', BirdParameters(i).Femaleintroductiontime, '; latency = ', num2str(LatencyToFirstSong(i)), ' mins or ', num2str(LatencyToFirstSong(i)*60), ' s']);
    end
end

LatencyToFirstSong(find((LatencyToFirstSong < 0) & (LatencyToFirstSong > -0.5))) = 0;
LatencyToFirstMotif(find((LatencyToFirstMotif < 0) & (LatencyToFirstMotif > -0.5))) = 0;

% Now for each bird, plot the latencies to first song for the different
% conditions - again do this separately for each mic type

% Get unique BirdNames to get an idea of how many birds there are
for i = 1:length(BirdParameters),
    BirdNames{i} = BirdParameters(i).BirdName;
    MicrophoneTypes{i} = BirdParameters(i).Microphone;
end

UniqueBirds = unique(BirdNames);

figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [318 66 1500 900]);
p = panel();
NumRows = ceil(length(UniqueBirds)/4)+2;
p.pack(NumRows, 4);
p.fontsize = 10;

Index = 0;
for i = 1:length(UniqueBirds),
    
    Matches = strmatch(UniqueBirds{i}, BirdNames);
    UniqueMicrophones = unique(MicrophoneTypes(Matches));
    
    for MicType = 1:length(UniqueMicrophones),
        % First initialize first song latency values to empty arrays for all
        % conditions (6 in all)
        Index = Index + 1;
        p(mod(Index-1,NumRows) + 1, ceil(Index/NumRows)).select();
        for j = 1:6,
            AcrossBirdFirstSongLatency{i}{MicType}{j} = [];
            AcrossBirdFirstMotifLatency{i}{MicType}{j} = [];
        end
        
        Matches = intersect(strmatch(UniqueBirds{i}, BirdNames), strmatch(UniqueMicrophones{MicType}, MicrophoneTypes));
        for j = 1:length(Matches),
            switch BirdParameters(Matches(j)).Condition
                case 'L0'
                    ConditionIndex = 1;
                case 'L1'
                    ConditionIndex = 2;
                case 'L2'
                    ConditionIndex = 3;
                case 'L3'
                    ConditionIndex = 4;
                case 'L4'
                    ConditionIndex = 5;
                case 'UN'
                    ConditionIndex = 6;
            end
            
            AcrossBirdFirstSongLatency{i}{MicType}{ConditionIndex}(end+1) = LatencyToFirstSong(Matches(j));
            AcrossBirdFirstMotifLatency{i}{MicType}{ConditionIndex}(end+1) = LatencyToFirstMotif(Matches(j));
        end
        
        title([UniqueBirds{i} ':' UniqueMicrophones{MicType}]);
        
        plot(cellfun(@median, AcrossBirdFirstMotifLatency{i}{MicType}), 'ro-');
        hold on;
        plot(cellfun(@median, AcrossBirdFirstSongLatency{i}{MicType}), 'bo-');
        for j = 1:length(AcrossBirdFirstMotifLatency{i}{MicType}),
            if (~isempty(AcrossBirdFirstMotifLatency{i}{MicType}{j}))
                plot(ones(1,2)*j, prctile(AcrossBirdFirstMotifLatency{i}{MicType}{j}, [25 75]), 'r');
                plot(ones(1,2)*j, prctile(AcrossBirdFirstSongLatency{i}{MicType}{j}, [25 75]), 'b');
            end
        end
        set(gca, 'XTick', 1:1:6, 'XTickLabel', {'L0' 'L1' 'L2' 'L3' 'L4' 'UN'}, 'FontSize', 10);
        
        if (ceil(Index/NumRows) == 1)
            ylabel('Latency (mins)');
        end
        if (mod(Index,NumRows) == 0)
            xlabel('Distance');
        end
        axis tight;
        Temp = axis;
        Temp = [0.5 6.5 0 1.04*Temp(4)];
        axis(Temp);
    end
end
OutputDir = '/home/raghav/StudentRelated/Harini';
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(OutputDir, 'LatencyToFirstSongsMotifs_SeparateMicrophones.png'), '-dpng');

figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [1487 250 450 850]);
p = panel();
p.pack(2, 1);
p.fontsize = 10;

for i = 1:length(UniqueBirds),
    Matches = strmatch(UniqueBirds{i}, BirdNames);
    UniqueMicrophones = unique(MicrophoneTypes(Matches));
    TempFirstMotifLatency = [];
    TempFirstSongLatency = [];
    for MicType = 1:length(UniqueMicrophones),
        TempFirstSongLatency(MicType, :) = cellfun(@median, AcrossBirdFirstSongLatency{i}{MicType});
        TempFirstMotifLatency(MicType, :) = cellfun(@median, AcrossBirdFirstMotifLatency{i}{MicType});
    end        
    
    if (size(TempFirstSongLatency,1) > 1)
        FirstSongLatency_AcrossBirds(i,:) = mean(TempFirstSongLatency);
        FirstMotifLatency_AcrossBirds(i,:) = mean(TempFirstMotifLatency);
    else
        FirstSongLatency_AcrossBirds(i,:) = (TempFirstSongLatency);
        FirstMotifLatency_AcrossBirds(i,:) = (TempFirstMotifLatency);
    end
end


for i = 1:size(FirstSongLatency_AcrossBirds,1),
    NonNaNVals = find(~isnan(FirstSongLatency_AcrossBirds(i,:)));
    p(1,1).select();
    hold on;
    plot(NonNaNVals, FirstSongLatency_AcrossBirds(i,NonNaNVals), 'ko-', 'Color', [0.8 0.8 0.8]);
    
    NonNaNVals = find(~isnan(FirstMotifLatency_AcrossBirds(i,:)));
    p(2,1).select();
    hold on;
    plot(NonNaNVals, FirstMotifLatency_AcrossBirds(i,NonNaNVals), 'ko-', 'Color', [0.8 0.8 0.8]);
end

p(1,1).select();
plot(nanmean(FirstSongLatency_AcrossBirds), 'ks-', 'MarkerSize', 8, 'LineWidth', 1.5);
for i = 1:size(FirstSongLatency_AcrossBirds, 2),
    errorbar(i,nanmean(FirstSongLatency_AcrossBirds(:,i)), nanstd(FirstSongLatency_AcrossBirds(:,i))/sqrt(length(find(~isnan(FirstSongLatency_AcrossBirds(:,i))))), 'ks', 'MarkerSize', 8, 'LineWidth', 1.5);
end
ylabel({'First song'; 'bout latency (mins)'});
% set(gca, 'YScale', 'log');
axis tight;
Temp = axis;
Temp = [0.5 6.5 0 1.02*Temp(4)];
axis(Temp);

p(2,1).select();
plot(nanmean(FirstMotifLatency_AcrossBirds), 'ks-', 'MarkerSize', 8, 'LineWidth', 1.5);
for i = 1:size(FirstMotifLatency_AcrossBirds, 2),
    errorbar(i,nanmean(FirstMotifLatency_AcrossBirds(:,i)), nanstd(FirstMotifLatency_AcrossBirds(:,i))/sqrt(length(find(~isnan(FirstMotifLatency_AcrossBirds(:,i))))), 'ks', 'MarkerSize', 8, 'LineWidth', 1.5);
end
ylabel({'First motif'; 'latency (mins)'});
set(gca, 'XTick', 1:1:6, 'XTickLabel', {'L0' 'L1' 'L2' 'L3' 'L4' 'UN'}, 'FontSize', 10);
xlabel('Distance');
% set(gca, 'YScale', 'log');
axis tight;
Temp = axis;
Temp = [0.5 6.5 0 1.02*Temp(4)];
axis(Temp);

p.fontsize = 14;
p.marginleft = 30;

OutputDir = '/home/raghav/StudentRelated/Harini';
set(gcf, 'PaperPositionMode', 'auto');
print(fullfile(OutputDir, 'LatencyToFirstSongsMotifs_GroupData.png'), '-dpng');

disp('Got bout nos');


