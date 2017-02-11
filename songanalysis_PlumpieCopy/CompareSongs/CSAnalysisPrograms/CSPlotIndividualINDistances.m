function [] = CSPlotIndividualINDistances(CSData)

FeatureAxisLabels = [{'MeanFrequency'} {'Mean frequency (Hz)'}; {'LogAmplitude'} {'Amplitude (dB)'}; {'Entropy'} {'Entropy'}; {'FundamentalFrequency'} {'Frequency (Hz)'}; {'MeanFrequency'} {'Mean frequency (Hz)'}];
% Algorithm
% For each day and for each bout, first get the first motif syllable. The
% syllables before this first motif syllable can be counted as Total
% Syllables before first motif syllable. In addition, the number of INs in
% these syllables can be counted for #of INs.

MotifSyllArray = cellstr(char(ones(length(CSData.AllLabels), 1)*double(CSData.MotifSyllLabels)));
INArray = cellstr(char(ones(length(CSData.AllLabels), 1)*double(CSData.INLabels)));
MotifInitiationSyllArray = cellstr(char(ones(length(CSData.AllLabels), 1)*double(CSData.MotifInitiationSyllLabels)));
    
% using cellfun so that i iterate over each element of the cell array.
% To use cellfun, all of the other inputs also have to be in the form
% of cell arrays of the same length - so the previous three lines
% convert file type, data dir and output dir - common parameters for
% all of the files into cell arrays
    
[INs, Motifs, Bouts] = cellfun(@CSIdentifyINs, CSData.AllLabels', MotifSyllArray, INArray, 'UniformOutput', 0);

% Now calculate mean distances of various INs to the last IN on that particular day for all days 
figure;

Colors = 'rgbcmk';
Symbols = 'o+d';

DistanceFeatures = [{'LogAmplitude'}; {'Duration'}; {'Entropy'}; {'MeanFrequency'}];

for i = 1:CSData.NoofDays,
    subplot(CSData.NoofDays, 1, i);
    hold on;
    
    for j = 1:length(DistanceFeatures),
        DistanceFeatureIndex(j) = strmatch(DistanceFeatures{j}, CSData.Data{i}.ToBeUsedFeatures, 'exact');
    end
    
    PosFromLast = sum(INs{i}.PosFromLast);
    
    LastINIndices = find(PosFromLast == -1);
    LastINFeats = CSData.AllFeats{i}(INs{i}.Indices(LastINIndices), DistanceFeatureIndex);
    
    CSData.DistanceFromLast{i} = [];
    for j = 1:length(INs{i}.Starts),
        Indices = INs{i}.Starts(j):INs{i}.Ends(j);
        Distances = pdist2(CSData.AllFeats{i}(Indices, DistanceFeatureIndex), mean(LastINFeats), 'mahalanobis', cov(LastINFeats));
        plot(-(length(Indices)):1:-1, Distances, 'ko-', 'MarkerSize', 2);
    end
end

for i = 1:CSData.NoofDays,
    subplot(CSData.NoofDays, 1, i);
    hold on;
    axis tight;
    Temp = axis;
    Temp(1) = Temp(1) - 0.5;
    Temp(2) = Temp(2) + 0.5;
    Temp(3) = 0;
    Temp(4) = Temp(4)*1.02;
    axis(Temp);
    plot(Temp(1:2), [2 2], 'b--', 'LineWidth', 2);
end

figure;
hold on;

Colors = 'rgbcmk';
Symbols = 'o+d';

DistanceFeatures = [{'LogAmplitude'}; {'Duration'}; {'Entropy'}; {'MeanFrequency'}];

for i = 1:CSData.NoofDays,
    
    for j = 1:length(DistanceFeatures),
        DistanceFeatureIndex(j) = strmatch(DistanceFeatures{j}, CSData.Data{i}.ToBeUsedFeatures, 'exact');
    end
    
    PosFromLast = sum(INs{i}.PosFromLast);
    
    LastINIndices = find(PosFromLast == -1);
    LastINFeats = CSData.AllFeats{i}(INs{i}.Indices(LastINIndices), DistanceFeatureIndex);
    
    CSData.DistanceFromLast{i} = [];
    for j = 1:length(INs{i}.Starts),
        Indices = INs{i}.Starts(j):INs{i}.Ends(j);
        Distances = pdist2(CSData.AllFeats{i}(Indices, DistanceFeatureIndex), mean(LastINFeats), 'mahalanobis', cov(LastINFeats));
        plot(-(length(Indices)):1:-1, Distances, [Colors(i), 'o-'], 'MarkerSize', 2);
    end
end

for i = 1:CSData.NoofDays,
    axis tight;
    Temp = axis;
    Temp(1) = Temp(1) - 0.5;
    Temp(2) = Temp(2) + 0.5;
    Temp(3) = 0;
    Temp(4) = Temp(4)*1.02;
    axis(Temp);
    plot(Temp(1:2), [2 2], 'k--', 'LineWidth', 2);
end

disp('Finished analysing data');