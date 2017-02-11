function [BoutLength] = MA_AnalyseSavedData_BoutLength(Parameters, TitleString)

%====== First analyse bout length for pre and post ========================
% Bout length is stored for each day in an array
% Row indicates a day and the first three columns are # of bouts, mean bout
% length and standard deviation of bout length for dir song. The next 3
% columns are the same for undir song. 
% First set of rows are for pre days and the next set of rows are for post
% days

OutputDir = '/home/raghav/HVC_MicrolesionDataFigures/PaperFigures/';

PrePostDays = [1 2; 2 3; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 2 5; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];

for ParameterNo = 1:length(Parameters),
    for i = 1:Parameters(ParameterNo).NoPreDays,
        DirBouts = cell2mat(Parameters(ParameterNo).PreDirBoutLens{i});
        UnDirBouts = cell2mat(Parameters(ParameterNo).PreUnDirBoutLens{i});

        BoutLength(ParameterNo).Dir.NumBouts(i) = length(DirBouts);
        BoutLength(ParameterNo).Dir.Mean(i) = mean(DirBouts);
        BoutLength(ParameterNo).Dir.STD(i) = std(DirBouts);
        BoutLength(ParameterNo).Dir.Median(i) = median(DirBouts);
        BoutLength(ParameterNo).Dir.IQR(i) = iqr(DirBouts);

        BoutLength(ParameterNo).UnDir.NumBouts(i) = length(UnDirBouts);
        BoutLength(ParameterNo).UnDir.Mean(i) = mean(UnDirBouts);
        BoutLength(ParameterNo).UnDir.STD(i) = std(UnDirBouts);
        BoutLength(ParameterNo).UnDir.Median(i) = median(UnDirBouts);
        BoutLength(ParameterNo).UnDir.IQR(i) = iqr(UnDirBouts);

        clear DirBouts UnDirBouts;
    end

    for i = 1:Parameters(ParameterNo).NoPostDays,
        DirBouts = cell2mat(Parameters(ParameterNo).PostDirBoutLens{i});
        UnDirBouts = cell2mat(Parameters(ParameterNo).PostUnDirBoutLens{i});

        BoutLength(ParameterNo).Dir.NumBouts(i + Parameters(ParameterNo).NoPreDays) = length(DirBouts);
        BoutLength(ParameterNo).Dir.Mean(i + Parameters(ParameterNo).NoPreDays) = mean(DirBouts);
        BoutLength(ParameterNo).Dir.STD(i + Parameters(ParameterNo).NoPreDays) = std(DirBouts);
        BoutLength(ParameterNo).Dir.Median(i + Parameters(ParameterNo).NoPreDays) = median(DirBouts);
        BoutLength(ParameterNo).Dir.IQR(i + Parameters(ParameterNo).NoPreDays) = iqr(DirBouts);

        BoutLength(ParameterNo).UnDir.NumBouts(i + Parameters(ParameterNo).NoPreDays) = length(UnDirBouts);
        BoutLength(ParameterNo).UnDir.Mean(i + Parameters(ParameterNo).NoPreDays) = mean(UnDirBouts);
        BoutLength(ParameterNo).UnDir.STD(i + Parameters(ParameterNo).NoPreDays) = std(UnDirBouts);
        BoutLength(ParameterNo).UnDir.Median(i + Parameters(ParameterNo).NoPreDays) = median(UnDirBouts);
        BoutLength(ParameterNo).UnDir.IQR(i + Parameters(ParameterNo).NoPreDays) = iqr(UnDirBouts);

        clear DirBouts UnDirBouts;
    end
end

%==========================================================================

for i = 1:length(BoutLength),
    DirMedians(i,1) = BoutLength(i).Dir.Median(PrePostDays(i,1));
    DirMedians(i,2) = BoutLength(i).Dir.Median(PrePostDays(i,2));
    
    UnDirMedians(i,1) = BoutLength(i).UnDir.Median(PrePostDays(i,1));
    UnDirMedians(i,2) = BoutLength(i).UnDir.Median(PrePostDays(i,2));
end

MA_PlotVsLesionSize([Parameters.PercentTotalHVCremaining], DirMedians, UnDirMedians, 'Ratio of Median Bout Length (Post/Pre)', OutputDir, 'SouthEast', 'MedianBoutLength', TitleString);

MA_PlotPrevsPost(DirMedians/1000, UnDirMedians/1000, 'Median Bout Length (s)', OutputDir, 'MedianBoutLength', TitleString);

MA_PlotDirVsUnDir([DirMedians(:,1) UnDirMedians(:,1)]/1000, [DirMedians(:,2) UnDirMedians(:,2)]/1000, 'Median Bout Length (s)', OutputDir, 'MedianBoutLength', TitleString);

disp('Finished plotting bout lengths');