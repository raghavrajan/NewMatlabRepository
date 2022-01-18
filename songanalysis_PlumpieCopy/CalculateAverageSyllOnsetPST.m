function [NormalizedSyllOnsetPST, Unit_Bird_SyllDetails] = CalculateAverageSyllOnsetPST(SyllOnsetCombinedData, SURecordingDetails, NeuronType, MinFR)

% For every unit and syllable, check if there are > 10 renditions of Dir
% and undir song and then calculate normalized smoothed pst - normalize to
% max of dir and undir

MinTrials = 10;

NormalizedSyllOnsetPST.Directed = [];
NormalizedSyllOnsetPST.UnDirected = [];

Index = 1;
for i = 1:length(SyllOnsetCombinedData),
    if (isempty(find(strfind(SURecordingDetails(SyllOnsetCombinedData{i}.DataSets(1)).Neurontype, NeuronType))))
        continue;
    end
    for j = 1:length(SyllOnsetCombinedData{i}.SmoothedPST),
        DirTrials = find(SyllOnsetCombinedData{i}.DirUnDir{j} == 1);
        UnDirTrials = find(SyllOnsetCombinedData{i}.DirUnDir{j} == 0);
        
        if ((length(DirTrials) >= MinTrials) & (length(UnDirTrials) >= MinTrials))
            DirSmoothedPST = mean(SyllOnsetCombinedData{i}.SmoothedPST{j}(DirTrials,:));
            UnDirSmoothedPST = mean(SyllOnsetCombinedData{i}.SmoothedPST{j}(UnDirTrials,:));
            
            NormalizedSyllOnsetPST.Directed(end+1,:) = DirSmoothedPST/(max([DirSmoothedPST(:); UnDirSmoothedPST(:)]));
            NormalizedSyllOnsetPST.UnDirected(end+1,:) = UnDirSmoothedPST/(max([DirSmoothedPST(:); UnDirSmoothedPST(:)]));
            Unit_Bird_SyllDetails.Unit_Syll(Index,:) = [i j];
            Unit_Bird_SyllDetails.Bird{Index} = SURecordingDetails(SyllOnsetCombinedData{i}.DataSets(1)).BirdName;
            Index = Index + 1;
        end
    end
end
            