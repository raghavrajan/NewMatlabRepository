function [NumSylls] = MA_AnalyseSavedData_NumSyllables(Parameters, TitleString, varargin)

if (nargin > 2)
    BirdIndices = varargin{1};
else
    BirdIndices = 1:1:length(Parameters);
end

%====== Analyse number of syllables ========================

OutputDir = '/home/raghav/HVC_MicrolesionDataFigures/PaperFigures/';
PrePostDays = [1 2; 2 3; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 2 5; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];

Parameters = Parameters(BirdIndices);
PrePostDays = PrePostDays(BirdIndices,:);

for ParameterNo = 1:length(Parameters),
    for i = 1:Parameters(ParameterNo).NoPreDays,
        BoutIndex = 1;
        for j = 1:length(Parameters(ParameterNo).PreDirBoutLens{i}),
            for k = 1:length(Parameters(ParameterNo).PreDirBoutLens{i}{j}),
                SyllableIndices = find((Parameters(ParameterNo).PreDirOnsets{i}{j} >= Parameters(ParameterNo).PreDirBoutOnsets{i}{j}(k)) & (Parameters(ParameterNo).PreDirOnsets{i}{j} <= Parameters(ParameterNo).PreDirBoutOffsets{i}{j}(k)));
                NumSylls(ParameterNo).Dir{i}(BoutIndex) = length(SyllableIndices);
                BoutIndex = BoutIndex + 1;
            end
        end

        BoutIndex = 1;
        for j = 1:length(Parameters(ParameterNo).PreUnDirBoutLens{i}),
            for k = 1:length(Parameters(ParameterNo).PreUnDirBoutLens{i}{j}),
                SyllableIndices = find((Parameters(ParameterNo).PreUnDirOnsets{i}{j} >= Parameters(ParameterNo).PreUnDirBoutOnsets{i}{j}(k)) & (Parameters(ParameterNo).PreUnDirOnsets{i}{j} <= Parameters(ParameterNo).PreUnDirBoutOffsets{i}{j}(k)));
                NumSylls(ParameterNo).UnDir{i}(BoutIndex) = length(SyllableIndices);
                BoutIndex = BoutIndex + 1;
            end
        end
    end


    for i = 1:Parameters(ParameterNo).NoPostDays,
        BoutIndex = 1;
        for j = 1:length(Parameters(ParameterNo).PostDirBoutLens{i}),
            for k = 1:length(Parameters(ParameterNo).PostDirBoutLens{i}{j}),
                SyllableIndices = find((Parameters(ParameterNo).PostDirOnsets{i}{j} >= Parameters(ParameterNo).PostDirBoutOnsets{i}{j}(k)) & (Parameters(ParameterNo).PostDirOnsets{i}{j} <= Parameters(ParameterNo).PostDirBoutOffsets{i}{j}(k)));
                NumSylls(ParameterNo).Dir{i + Parameters(ParameterNo).NoPreDays}(BoutIndex) = length(SyllableIndices);
                BoutIndex = BoutIndex + 1;
            end
        end

        BoutIndex = 1;
        for j = 1:length(Parameters(ParameterNo).PostUnDirBoutLens{i}),
            for k = 1:length(Parameters(ParameterNo).PostUnDirBoutLens{i}{j}),
                SyllableIndices = find((Parameters(ParameterNo).PostUnDirOnsets{i}{j} >= Parameters(ParameterNo).PostUnDirBoutOnsets{i}{j}(k)) & (Parameters(ParameterNo).PostUnDirOnsets{i}{j} <= Parameters(ParameterNo).PostUnDirBoutOffsets{i}{j}(k)));
                NumSylls(ParameterNo).UnDir{i + Parameters(ParameterNo).NoPreDays}(BoutIndex) = length(SyllableIndices);
                BoutIndex = BoutIndex + 1;
            end
        end
    end
end

for i = 1:length(NumSylls),
    DirMedians(i,1) = median(NumSylls(i).Dir{PrePostDays(i,1)});
    DirMedians(i,2) = median(NumSylls(i).Dir{PrePostDays(i,2)});
    
    UnDirMedians(i,1) = median(NumSylls(i).UnDir{PrePostDays(i,1)});
    UnDirMedians(i,2) = median(NumSylls(i).UnDir{PrePostDays(i,2)});
end

MA_PlotVsLesionSize([Parameters.PercentTotalHVCremaining], DirMedians, UnDirMedians, 'Median # of syllables per bout (Post/Pre)', OutputDir, 'SouthEast', 'NumSyllsPerBout', TitleString);

MA_PlotPrevsPost(DirMedians, UnDirMedians, 'Median # of syllables per bout', OutputDir, 'NumSyllsPerBout', TitleString);

MA_PlotDirVsUnDir([DirMedians(:,1) UnDirMedians(:,1)], [DirMedians(:,2) UnDirMedians(:,2)], 'Median # of syllables per bout', OutputDir, 'NumSyllsPerBout', TitleString);

disp('Finished plotting num syllables per bout');