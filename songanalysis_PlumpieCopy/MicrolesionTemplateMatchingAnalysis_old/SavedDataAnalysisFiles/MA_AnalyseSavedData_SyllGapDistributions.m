function [SyllGapStats] = MA_AnalyseSavedData_SyllGapDistributions(Parameters, TitleString, varargin)

if (nargin > 2)
    BirdIndices = varargin{1};
else
    BirdIndices = 1:1:length(Parameters);
end

OutputDir = '/home/raghav/HVC_MicrolesionDataFigures/PaperFigures/';

PrePostDays = [1 2; 2 3; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 2 5; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];

Parameters = Parameters(BirdIndices);
PrePostDays = PrePostDays(BirdIndices, :);


Params.NoofBins = 50;
Params.DurEdges = linspace(-2.5, -0.3, Params.NoofBins);
min_int = 7;
min_dur = 7;
                    
for ParameterNo = 1:length(Parameters),
    fprintf('%s >>', Parameters(ParameterNo).BirdName);
    NoteFileDir = [OutputDir, filesep, Parameters(ParameterNo).BirdName, '.NoteFiles', filesep];
    if (~exist(NoteFileDir, 'dir'))
        mkdir(NoteFileDir);
    end
    for i = 1:Parameters(ParameterNo).NoPreDays,
        fprintf(' Pre Day %i >>', i);
        fprintf(' Dir >>');
        
        DirSyllDurs = [];
        DirGapDurs = [];
        BoutIndex = 1;
        for j = 1:length(Parameters(ParameterNo).PreDirBoutLens{i}),
            [RawData, Fs] = ASSLGetRawData(Parameters(ParameterNo).PreDataDir{i}, Parameters(ParameterNo).PreDirSongFileNames{i}{j}, Parameters(ParameterNo).FileType, 0);
            for k = 1:length(Parameters(ParameterNo).PreDirBoutLens{i}{j}),
                SongBout = RawData(ceil(Parameters(ParameterNo).PreDirBoutOnsets{i}{j}(k)*Fs/1000):floor(Parameters(ParameterNo).PreDirBoutOffsets{i}{j}(k)*Fs/1000));
                
                NoteFile = [NoteFileDir, Parameters(ParameterNo).PreDirSongFileNames{i}{j}, '.Bout', num2str(k), '.not.mat'];
                if (exist(NoteFile, 'file'))
                    Temp = load(NoteFile);
                    Onsets = Temp.onsets;
                    Offsets = Temp.offsets;
                else
                    BoutOnset = Parameters(ParameterNo).PreDirBoutOnsets{i}{j}(k);
                    BoutOffset = Parameters(ParameterNo).PreDirBoutOffsets{i}{j}(k);
                    [Onsets, Offsets] = MA_AnalyseSavedData_SegmentFilesAronovFee(SongBout, Fs, min_int, min_dur, NoteFile, BoutOnset, BoutOffset);
                end
            
                SyllDurs = Offsets - Onsets;
                GapDurs = Onsets(2:end) - Offsets(1:end-1);
                
                DirSyllDurs = [DirSyllDurs; SyllDurs(:)];
                DirGapDurs = [DirGapDurs; GapDurs(:)];

                clear SyllDurs GapDurs;
                BoutIndex = BoutIndex + 1;
            end
        end

        fprintf(' Undir >>');
        
        UnDirSyllDurs = [];
        UnDirGapDurs = [];
        BoutIndex = 1;
        for j = 1:length(Parameters(ParameterNo).PreUnDirBoutLens{i}),
            [RawData, Fs] = ASSLGetRawData(Parameters(ParameterNo).PreDataDir{i}, Parameters(ParameterNo).PreUnDirSongFileNames{i}{j}, Parameters(ParameterNo).FileType, 0);
            for k = 1:length(Parameters(ParameterNo).PreUnDirBoutLens{i}{j}),
                SongBout = RawData(ceil(Parameters(ParameterNo).PreUnDirBoutOnsets{i}{j}(k)*Fs/1000):floor(Parameters(ParameterNo).PreUnDirBoutOffsets{i}{j}(k)*Fs/1000));
                
                NoteFile = [NoteFileDir, Parameters(ParameterNo).PreUnDirSongFileNames{i}{j}, '.Bout', num2str(k), '.not.mat'];
                if (exist(NoteFile, 'file'))
                    Temp = load(NoteFile);
                    Onsets = Temp.onsets;
                    Offsets = Temp.offsets;
                else
                    BoutOnset = Parameters(ParameterNo).PreUnDirBoutOnsets{i}{j}(k);
                    BoutOffset = Parameters(ParameterNo).PreUnDirBoutOffsets{i}{j}(k);
                    [Onsets, Offsets] = MA_AnalyseSavedData_SegmentFilesAronovFee(SongBout, Fs, min_int, min_dur, NoteFile, BoutOnset, BoutOffset);
                end

                SyllDurs = Offsets - Onsets;
                GapDurs = Onsets(2:end) - Offsets(1:end-1);
                
                UnDirSyllDurs = [UnDirSyllDurs; SyllDurs(:)];
                UnDirGapDurs = [UnDirGapDurs; GapDurs(:)];

                clear SyllDurs GapDurs;
                BoutIndex = BoutIndex + 1;
            end
        end

        SyllGapStats(ParameterNo).DirSyllDurs{i} = DirSyllDurs;
        SyllGapStats(ParameterNo).DirGapDurs{i} = DirGapDurs;
        SyllGapStats(ParameterNo).UnDirSyllDurs{i} = UnDirSyllDurs;
        SyllGapStats(ParameterNo).UnDirGapDurs{i} = UnDirGapDurs;
        
        DirSyllDurs = log10(DirSyllDurs/1000);
        DirGapDurs = log10(DirGapDurs/1000);
        UnDirSyllDurs = log10(UnDirSyllDurs/1000);
        UnDirGapDurs = log10(UnDirGapDurs/1000);

        SyllGapStats(ParameterNo).DirSyllDurHist(i,:) = histc(DirSyllDurs, Params.DurEdges)/length(DirSyllDurs);
        SyllGapStats(ParameterNo).UnDirSyllDurHist(i,:) = histc(UnDirSyllDurs, Params.DurEdges)/length(UnDirSyllDurs);
        SyllGapStats(ParameterNo).DirGapDurHist(i,:) = histc(DirGapDurs, Params.DurEdges)/length(DirGapDurs);
        SyllGapStats(ParameterNo).UnDirGapDurHist(i,:) = histc(UnDirGapDurs, Params.DurEdges)/length(UnDirGapDurs);    
    end


    for i = 1:Parameters(ParameterNo).NoPostDays,
        fprintf(' Post Day %i >>', i);
        fprintf(' Dir >>');
        
        DirSyllDurs = [];
        DirGapDurs = [];
        BoutIndex = 1;
        for j = 1:length(Parameters(ParameterNo).PostDirBoutLens{i}),
            [RawData, Fs] = ASSLGetRawData(Parameters(ParameterNo).PostDataDir{i}, Parameters(ParameterNo).PostDirSongFileNames{i}{j}, Parameters(ParameterNo).FileType, 0);
            for k = 1:length(Parameters(ParameterNo).PostDirBoutLens{i}{j}),
                SongBout = RawData(ceil(Parameters(ParameterNo).PostDirBoutOnsets{i}{j}(k)*Fs/1000):floor(Parameters(ParameterNo).PostDirBoutOffsets{i}{j}(k)*Fs/1000));
                
                NoteFile = [NoteFileDir, Parameters(ParameterNo).PostDirSongFileNames{i}{j}, '.Bout', num2str(k), '.not.mat'];
                if (exist(NoteFile, 'file'))
                    Temp = load(NoteFile);
                    Onsets = Temp.onsets;
                    Offsets = Temp.offsets;
                else
                    BoutOnset = Parameters(ParameterNo).PostDirBoutOnsets{i}{j}(k);
                    BoutOffset = Parameters(ParameterNo).PostDirBoutOffsets{i}{j}(k);
                    [Onsets, Offsets] = MA_AnalyseSavedData_SegmentFilesAronovFee(SongBout, Fs, min_int, min_dur, NoteFile, BoutOnset, BoutOffset);
                end

                SyllDurs = Offsets - Onsets;
                GapDurs = Onsets(2:end) - Offsets(1:end-1);
                
                DirSyllDurs = [DirSyllDurs; SyllDurs(:)];
                DirGapDurs = [DirGapDurs; GapDurs(:)];

                clear SyllDurs GapDurs;
                BoutIndex = BoutIndex + 1;
            end
        end

        fprintf(' Undir\n');
        
        UnDirSyllDurs = [];
        UnDirGapDurs = [];
        BoutIndex = 1;
        for j = 1:length(Parameters(ParameterNo).PostUnDirBoutLens{i}),
            [RawData, Fs] = ASSLGetRawData(Parameters(ParameterNo).PostDataDir{i}, Parameters(ParameterNo).PostUnDirSongFileNames{i}{j}, Parameters(ParameterNo).FileType, 0);
            for k = 1:length(Parameters(ParameterNo).PostUnDirBoutLens{i}{j}),
                SongBout = RawData(ceil(Parameters(ParameterNo).PostUnDirBoutOnsets{i}{j}(k)*Fs/1000):floor(Parameters(ParameterNo).PostUnDirBoutOffsets{i}{j}(k)*Fs/1000));
                
                NoteFile = [NoteFileDir, Parameters(ParameterNo).PostUnDirSongFileNames{i}{j}, '.Bout', num2str(k), '.not.mat'];
                if (exist(NoteFile, 'file'))
                    Temp = load(NoteFile);
                    Onsets = Temp.onsets;
                    Offsets = Temp.offsets;
                else
                    BoutOnset = Parameters(ParameterNo).PostUnDirBoutOnsets{i}{j}(k);
                    BoutOffset = Parameters(ParameterNo).PostUnDirBoutOffsets{i}{j}(k);
                    [Onsets, Offsets] = MA_AnalyseSavedData_SegmentFilesAronovFee(SongBout, Fs, min_int, min_dur, NoteFile, BoutOnset, BoutOffset);
                end

                SyllDurs = Offsets - Onsets;
                GapDurs = Onsets(2:end) - Offsets(1:end-1);
                
                UnDirSyllDurs = [UnDirSyllDurs; SyllDurs(:)];
                UnDirGapDurs = [UnDirGapDurs; GapDurs(:)];
             
                clear SyllDurs GapDurs;
                BoutIndex = BoutIndex + 1;
            end
        end

        SyllGapStats(ParameterNo).DirSyllDurs{i + Parameters(ParameterNo).NoPreDays} = DirSyllDurs;
        SyllGapStats(ParameterNo).DirGapDurs{i + Parameters(ParameterNo).NoPreDays} = DirGapDurs;
        SyllGapStats(ParameterNo).UnDirSyllDurs{i + Parameters(ParameterNo).NoPreDays} = UnDirSyllDurs;
        SyllGapStats(ParameterNo).UnDirGapDurs{i + Parameters(ParameterNo).NoPreDays} = UnDirGapDurs;

        DirSyllDurs = log10(DirSyllDurs/1000);
        DirGapDurs = log10(DirGapDurs/1000);
        UnDirSyllDurs = log10(UnDirSyllDurs/1000);
        UnDirGapDurs = log10(UnDirGapDurs/1000);

        SyllGapStats(ParameterNo).DirSyllDurHist(i + Parameters(ParameterNo).NoPreDays,:) = histc(DirSyllDurs, Params.DurEdges)/length(DirSyllDurs);
        SyllGapStats(ParameterNo).UnDirSyllDurHist(i + Parameters(ParameterNo).NoPreDays,:) = histc(UnDirSyllDurs, Params.DurEdges)/length(UnDirSyllDurs);
        SyllGapStats(ParameterNo).DirGapDurHist(i + Parameters(ParameterNo).NoPreDays,:) = histc(DirGapDurs, Params.DurEdges)/length(DirGapDurs);
        SyllGapStats(ParameterNo).UnDirGapDurHist(i + Parameters(ParameterNo).NoPreDays,:) = histc(UnDirGapDurs, Params.DurEdges)/length(UnDirGapDurs);    
    end

    % Now to store the entropy-based measure of variability (Goldberg and Fee,
    % JNeurophys, 2011) in a matrix - each row corresponds to a day and there
    % are two values in each row, the first one corresponding to dir and the
    % second one to undir. Same for gaps.

    for i = 1:size(SyllGapStats(ParameterNo).DirSyllDurHist,1),
        Temp = SyllGapStats(ParameterNo).DirSyllDurHist(i,:).*log10(SyllGapStats(ParameterNo).DirSyllDurHist(i,:));
        SyllGapStats(ParameterNo).SyllVar.Dir(i) = -sum(Temp(~isnan(Temp)))/log10(Params.NoofBins);

        Temp = SyllGapStats(ParameterNo).UnDirSyllDurHist(i,:).*log10(SyllGapStats(ParameterNo).UnDirSyllDurHist(i,:));
        SyllGapStats(ParameterNo).SyllVar.UnDir(i) = -sum(Temp(~isnan(Temp)))/log10(Params.NoofBins);

        Temp = SyllGapStats(ParameterNo).DirGapDurHist(i,:).*log10(SyllGapStats(ParameterNo).DirGapDurHist(i,:));
        SyllGapStats(ParameterNo).GapVar.Dir(i) = -sum(Temp(~isnan(Temp)))/log10(Params.NoofBins);

        Temp = SyllGapStats(ParameterNo).UnDirGapDurHist(i,:).*log10(SyllGapStats(ParameterNo).UnDirGapDurHist(i,:));
        SyllGapStats(ParameterNo).GapVar.UnDir(i) = -sum(Temp(~isnan(Temp)))/log10(Params.NoofBins);
    end
end

%==========================================================================
% First plot the syllable entropy changes like Goldberg and
% Fee, J.Neurophys. 2011

clear DirMedians UnDirMedians;
for i = 1:length(SyllGapStats),
    DirMedians(i,1) = median(SyllGapStats(i).SyllVar.Dir(PrePostDays(i,1)));
    DirMedians(i,2) = median(SyllGapStats(i).SyllVar.Dir(PrePostDays(i,2)));
    
    UnDirMedians(i,1) = median(SyllGapStats(i).SyllVar.UnDir(PrePostDays(i,1)));
    UnDirMedians(i,2) = median(SyllGapStats(i).SyllVar.UnDir(PrePostDays(i,2)));
end

MA_PlotVsLesionSize([Parameters.PercentTotalHVCremaining], DirMedians, UnDirMedians, 'Syllable variability (Post/Pre)', OutputDir, 'NorthWest', 'SyllDistVar_Entropy', TitleString);

MA_PlotPreVsPost(DirMedians, UnDirMedians, 'Syllable variability', OutputDir, 'SyllDistVar_Entropy', TitleString);

MA_PlotDirVsUnDir([DirMedians(:,1) UnDirMedians(:,1)], [DirMedians(:,2) UnDirMedians(:,2)], 'Syllable variability', OutputDir, 'SyllDistVar_Entropy', TitleString);
%==========================================================================
% Next plot the gap entropy changes like Goldberg and
% Fee, J.Neurophys. 2011

clear DirMedians UnDirMedians;

for i = 1:length(SyllGapStats),
    DirMedians(i,1) = median(SyllGapStats(i).GapVar.Dir(PrePostDays(i,1)));
    DirMedians(i,2) = median(SyllGapStats(i).GapVar.Dir(PrePostDays(i,2)));
    
    UnDirMedians(i,1) = median(SyllGapStats(i).GapVar.UnDir(PrePostDays(i,1)));
    UnDirMedians(i,2) = median(SyllGapStats(i).GapVar.UnDir(PrePostDays(i,2)));
end

MA_PlotVsLesionSize([Parameters.PercentTotalHVCremaining], DirMedians, UnDirMedians, 'Gap variability (Post/Pre)', OutputDir, 'NorthWest', 'GapDistVar_Entropy', TitleString);

MA_PlotPreVsPost(DirMedians, UnDirMedians, 'Gap variability', OutputDir, 'GapDistVar_Entropy', TitleString);

MA_PlotDirVsUnDir([DirMedians(:,1) UnDirMedians(:,1)], [DirMedians(:,2) UnDirMedians(:,2)], 'Gap variability', OutputDir, 'GapDistVar_Entropy', TitleString);
%==========================================================================
% Next find number of peaks in the syllable duration histogram like Naie
% and Hahnloser, J.Neurophys. 2011

clear DirMedians UnDirMedians;
BinSize = 5; % ms
Edges = 0:BinSize:500;

SmoothWins = (10:5:50); % in ms
Thresholds = 0.001:0.001:0.05; 

Fid = fopen([OutputDir, 'NumSyllablePeaks.SensititivityAnalysis.txt'], 'w');
% Print column labels
fprintf(Fid, '#\tSmoothing Window size (ms)\tThreshold for peak detection\tMean # of peaks (Pre)\t \tSEM # of peaks (Pre)\t \tMean # of peaks (Post)\t \tSEM # of peaks (Post)\t \tCorrelation between Post/Pre ratio and lesion size\t \t \t \t');
fprintf(Fid, 'Difference between # of peaks - Dir - Pre vs. Post (Sign rank test p)\tDifference between # of peaks - UnDir - Pre vs. Post (Sign rank test p)\t');
fprintf(Fid, 'Difference between # of peaks - Pre - Dir vs. UnDir (Sign rank test p)\tDifference between # of peaks - Post - Dir vs. UnDir (Sign rank test p)\t');
fprintf(Fid, '\n');

fprintf(Fid, '#\t \t \tDir\tUnDir\tDir\tUnDir\tDir\tUnDir\tDir\tUnDir\tDir (r)\tDir (p)\tUnDir (r)\tUnDir (p)\t');
fprintf(Fid, ' \t \t');
fprintf(Fid, ' \t \t');
fprintf(Fid, '\n');

for i = 1:length(SmoothWins), % ms
    for j = 1:length(Thresholds),

        SmoothWinSize = SmoothWins(i);
        Thresh = Thresholds(j);
        
        if ((SmoothWinSize == 15) && ((Thresh - 0.01) < 0.0000001))
            PlotYesNo = 1;
        else
            PlotYesNo = 0;
        end

        Smooth = ones(round(SmoothWinSize/BinSize), 1)/(SmoothWinSize/BinSize);

        clear SyllDurHist;

        for k = 1:length(SyllGapStats),
            SyllDurHist(1,:) = conv(histc(SyllGapStats(k).DirSyllDurs{PrePostDays(k,1)}, Edges), Smooth, 'same');
            SyllDurHist(2,:) = conv(histc(SyllGapStats(k).DirSyllDurs{PrePostDays(k,2)}, Edges), Smooth, 'same');
            SyllDurHist(3,:) = conv(histc(SyllGapStats(k).UnDirSyllDurs{PrePostDays(k,1)}, Edges), Smooth, 'same');
            SyllDurHist(4,:) = conv(histc(SyllGapStats(k).UnDirSyllDurs{PrePostDays(k,2)}, Edges), Smooth, 'same');

            for RowNo = 1:size(SyllDurHist, 1),
                SyllDurHist(RowNo,:) = SyllDurHist(RowNo,:)/sum(SyllDurHist(RowNo,:));
            end

%             DirMedians(i, 1) = length(findpeaks(abs(diff(SyllDurHist(1,:))), 'MINPEAKHEIGHT', Thresh));
%             DirMedians(i, 2) = length(findpeaks(abs(diff(SyllDurHist(2,:))), 'MINPEAKHEIGHT', Thresh));
%             UnDirMedians(i, 1) = length(findpeaks(abs(diff(SyllDurHist(3,:))), 'MINPEAKHEIGHT', Thresh));
%             UnDirMedians(i, 2) = length(findpeaks(abs(diff(SyllDurHist(4,:))), 'MINPEAKHEIGHT', Thresh));
            DirMedians(k, 1) = length(findpeaks((SyllDurHist(1,:)), 'MINPEAKHEIGHT', Thresh));
            DirMedians(k, 2) = length(findpeaks((SyllDurHist(2,:)), 'MINPEAKHEIGHT', Thresh));
            UnDirMedians(k, 1) = length(findpeaks((SyllDurHist(3,:)), 'MINPEAKHEIGHT', Thresh));
            UnDirMedians(k, 2) = length(findpeaks((SyllDurHist(4,:)), 'MINPEAKHEIGHT', Thresh));

        end

        MeanNumPeaks_Pre_Dir(i, j) = mean(DirMedians(:,1));
        MeanNumPeaks_Post_Dir(i, j) = mean(DirMedians(:,2));
        MeanNumPeaks_Pre_UnDir(i, j) = mean(UnDirMedians(:,1));
        MeanNumPeaks_Post_UnDir(i, j) = mean(UnDirMedians(:,2));
        
        SEMNumPeaks_Pre_Dir(i, j) = std(DirMedians(:,1))/sqrt(size(DirMedians,1));
        SEMNumPeaks_Post_Dir(i, j) = std(DirMedians(:,2))/sqrt(size(DirMedians,1));
        SEMNumPeaks_Pre_UnDir(i, j) = std(UnDirMedians(:,1))/sqrt(size(DirMedians,1));
        SEMNumPeaks_Post_UnDir(i, j) = std(UnDirMedians(:,2))/sqrt(size(DirMedians,1));
        
        
        [Dir_r, Dir_p, UnDir_r, UnDir_p] = MA_PlotVsLesionSize([Parameters.PercentTotalHVCremaining], DirMedians, UnDirMedians, '# of peaks in syll duration histogram (Post/Pre)', OutputDir, 'NorthWest', ['SmoothWinSize', num2str(SmoothWinSize), 'ms.NumPeaks_SyllDist'], TitleString, PlotYesNo);
        
        Dir_LesionSize_Corr_r(i, j) = Dir_r(1,2);
        Dir_LesionSize_Corr_p(i, j) = Dir_p(1,2);
        
        UnDir_LesionSize_Corr_r(i, j) = UnDir_r(1,2);
        UnDir_LesionSize_Corr_p(i, j) = UnDir_p(1,2);
        
        [PrePost_Dir_p(i, j), PrePost_UnDir_p(i, j)] = MA_PlotPreVsPost(DirMedians, UnDirMedians, '# of peaks in syll dur histogram', OutputDir, ['SmoothWinSize', num2str(SmoothWinSize), 'ms.NumPeaks_SyllDist'], TitleString, PlotYesNo);

        [DirUnDir_Pre_p(i, j), DirUnDir_Post_p(i, j)] = MA_PlotDirVsUnDir([DirMedians(:,1) UnDirMedians(:,1)], [DirMedians(:,2) UnDirMedians(:,2)], '# of peaks in syll dur histogram', OutputDir, ['SmoothWinSize', num2str(SmoothWinSize), 'ms.NumPeaks_SyllDist'], TitleString, PlotYesNo);

        fprintf(Fid, '%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t', i, SmoothWinSize, Thresh, MeanNumPeaks_Pre_Dir(i, j), MeanNumPeaks_Pre_UnDir(i, j), SEMNumPeaks_Pre_Dir(i, j), SEMNumPeaks_Pre_UnDir(i, j), MeanNumPeaks_Post_Dir(i, j), MeanNumPeaks_Post_UnDir(i, j), SEMNumPeaks_Post_Dir(i, j), SEMNumPeaks_Post_UnDir(i, j), Dir_LesionSize_Corr_r(i, j), Dir_LesionSize_Corr_p(i, j), UnDir_LesionSize_Corr_r(i, j), UnDir_LesionSize_Corr_p(i, j));
        fprintf(Fid, '%f\t%f', PrePost_Dir_p(i, j), PrePost_UnDir_p(i, j));
        fprintf(Fid, '%f\t%f', DirUnDir_Pre_p(i, j), DirUnDir_Post_p(i, j));
        fprintf(Fid, '\n');
    end
end

fclose(Fid);
disp('Finished analysing syllable gap statistics');