function [] = AutomaticallyIdentifyingMotifs(InterBoutInterval)

% File name of the file used to get the data from
SongDetailsFile = '/data/raghav/AutomaticMotifIdentification/ContinuousDataSongAnalysis_BirdDetails_ForAutomaticMotifIdentificationRelatedAnalysis.csv';
OutputDir = '/data/raghav/AutomaticMotifIdentification/SavedDataDir';
SavedDataDir = '/data/raghav/AutomaticMotifIdentification/SavedDataDir';

%% Get all data
[BirdParameters, Flag] = ProcessSongData_IntoBouts(SongDetailsFile, InterBoutInterval, SavedDataDir);

for i = 1:length(BirdParameters),
    [AllLabels{i}] = CombineContinuousDataNoteFiles(BirdParameters(i).DataDirectory, BirdParameters(i).SongFileNames, fullfile(BirdParameters(i).DataDirectory, 'ASSLNoteFiles'), BirdParameters(i).FileType);
    
    if (isfield(BirdParameters(i), 'Femalecalllabel'))
        FemaleINOUTEvents = regexp(AllLabels{i}, ['[', BirdParameters(i).Femalecalllabel, BirdParameters(i).FemaleinLabel, BirdParameters(i).FemaleoutLabel, BirdParameters(i).Closedoorlabel, ']']);
        if (~isempty(FemaleINOUTEvents))
            AllLabels{i}(FemaleINOUTEvents) = [];
        end
    end
    
    AllLabels{i} = AllLabels{i}(:);
end

disp('Finished');