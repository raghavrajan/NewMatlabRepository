function [Output] = Prasanth_AnalyzeNormalSongData(BirdDetailsTextFile, PlotOption, MinTrialNo, InterboutInterval, varargin)

% ========== Various options for plotting and help ========================

if (strfind(PlotOption, 'help'))
    disp('Usage: Prasanth_AnalyzeNormalSongData(CSVTextFile, PlotOption, MinTrialNo);');
    disp('PlotOption string can be any of the following');
end

% First process data from CSV file
[BirdParameters] = ProcessSongData(BirdDetailsTextFile, InterboutInterval);

% One more pre-processing step. I currently have a list of bouts and I also
% have a separate list of syllables that is compiled across all files. I
% have to have a set of indices that allow me to know which bout each
% syllable in the list of syllables belongs to
for i = 1:length(BirdParameters),
    BirdParameters(i).SyllableListBoutNum = zeros(size(BirdParameters(i).SyllableData,1), 1);
    TotalSyllNo = 0;
    for j = 1:length(BirdParameters(i).SongFileNames),
        BoutIndices = find(BirdParameters(i).Bouts(:,3) == j);
        for k = 1:length(BoutIndices),
            BirdParameters(i).SyllableListBoutNum((TotalSyllNo + BirdParameters(i).Bouts(BoutIndices(k), 1)):(TotalSyllNo + BirdParameters(i).Bouts(BoutIndices(k), 2))) = BoutIndices(k);
        end
        TotalSyllNo = TotalSyllNo + length(BirdParameters(i).NoteInfo{j}.labels);
    end
end

% Now Prasanth has two kinds of note files - one for the position and
% syllable feature analyis and one for the bout labels which can be used to
% analyze the properties of intervals
% These two have to be loaded up separately now.

disp('Loading up extra note data ...');
for i = 1:length(BirdParameters),
    BoutNotes_FileName = [BirdParameters(i).BirdName, '.', BirdParameters(i).DataLabel, '.', BirdParameters(i).Condition, '.', BirdParameters(i).Microphone, '.BoutNotes.mat'];
    if (exist(BoutNotes_FileName, 'file'))
        load(BoutNotes_FileName);
        BirdParameters(i).BoutNoteInfo = BoutNoteInfo;
    else
        [BirdParameters(i).BoutNoteInfo] = LSINA_LoadNoteFileInfo_ForPrasanth(BirdParameters(i));
        BoutNoteInfo = BirdParameters(i).BoutNoteInfo;
        save(BoutNotes_FileName, 'BoutNoteInfo');
    end
end

% ========= Analysis portion ==============================================

if (~isempty(which(['Prasanth_', PlotOption, '.m'])))
    if (nargin <= 4)
        Output = eval(['Prasanth_', PlotOption, '(BirdParameters, MinTrialNo);']);
    else
        Inputs = varargin{1};
        Output = eval(['Prasanth_', PlotOption, '(BirdParameters, MinTrialNo, Inputs);']);
    end
else
    disp(['Script Prasanth_', PlotOption, '.m does not exist']);
    return;
end    


disp('Finished');
