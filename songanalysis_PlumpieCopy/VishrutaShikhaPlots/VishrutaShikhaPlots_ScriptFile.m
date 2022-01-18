function [] = VishrutaShikhaPlots_ScriptFile(CSVTextFile, PlotOption)

[BirdParameters, Flag] = ProcessSongData(CSVTextFile, 2000); % 2000 here refers to inter-bout interval
if (sum(Flag) > 0)
    return;
end

Colours = 'rgbcmk';
Symbols = '+o<sd';

% One more pre-processing step. I currently have a list of bouts and I also
% have a separate list of syllables that is compiled across all files. I
% have to have a set of indices that allow me to know which bout each
% syllable in the list of syllables belongs to
for i = 1:length(BirdParameters),
    BirdParameters(i).SyllableListBoutNum = zeros(size(BirdParameters(i).SyllableData,1), 1);
    if (BirdParameters(i).Continuousdata == 0)
        TotalSyllNo = 0;
        for j = 1:length(BirdParameters(i).SongFileNames),
            BoutIndices = find(BirdParameters(i).Bouts(:,3) == j);
            FileSyllNo = 0;
            for k = 1:length(BoutIndices),
                BirdParameters(i).SyllableListBoutNum((TotalSyllNo + BirdParameters(i).Bouts(BoutIndices(k), 1)):(TotalSyllNo + BirdParameters(i).Bouts(BoutIndices(k), 2))) = BoutIndices(k);
                FileSyllNo = FileSyllNo + (BirdParameters(i).Bouts(BoutIndices(k), 2) - BirdParameters(i).Bouts(BoutIndices(k), 1) + 1);
            end
            TotalSyllNo = TotalSyllNo + FileSyllNo;
        end
    else
        for j = 1:size(BirdParameters(i).Bouts,1),
            BirdParameters(i).SyllableListBoutNum(BirdParameters(i).Bouts(j,1):BirdParameters(i).Bouts(j,2)) = j;
        end
    end
end

% Made one common way of plotting the results of any analysis. Basically,
% there is a script file called Harini_<<PlotOption>>.m and that is called
% to do the appropriate analysis. The first few lines are to check if the
% file exists - if not, then program returns with a message that such a
% program does not exist. 
% All scripts also take the same inputs to make things uniform

if (~isempty(which(['VishrutaShikha_', PlotOption, '.m'])))
    eval(['VishrutaShikha_', PlotOption, '(BirdParameters);']);
else
    disp(['Script VishrutaShikha_', PlotOption, '.m does not exist']);
    return;
end

disp('Finished plotting');