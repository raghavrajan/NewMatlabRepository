function [] = PlotSyllTransProbabilitiesAutoLabel(Dir, Dates, Motif, OutlierThreshold)

PresentDir = pwd;
cd(Dir);

for i = 1:length(Dates), 
    for j = 1:length(Motif),
        AnalysisFiles = dir(['*', Dates{i}, '*.motifseqmatches.', Motif(j), '.*.mat']);
        AnalysedData{i}{j} = load(AnalysisFiles(1).name);
        TempData = [];
        for k = 1:length(AnalysedData{i}{j}.UnDirBout.MaxBoutSeqMatch),
            TempData = [TempData findpeaks(AnalysedData{i}{j}.UnDirBout.MaxBoutSeqMatch{k}, 'MINPEAKHEIGHT', 0)];
        end
        for k = 1:length(AnalysedData{i}{j}.DirBout.MaxBoutSeqMatch),
            TempData = [TempData findpeaks(AnalysedData{i}{j}.DirBout.MaxBoutSeqMatch{k}, 'MINPEAKHEIGHT', 0)];
        end
        TempData = sort(TempData);
        MeanAnalysedData{i}{j} = mean(TempData);
        MedianAnalysedData{i}{j} = median(TempData);
        STDAnalysedData{i}{j} = std(TempData);
        MADAnalysedData{i}{j} = mad(TempData);
        Median_MADAnalysedData{i}{j} = mad(TempData, 1);
        Percentile_90_AnalysedData{i}{j} = TempData(round(0.9*length(TempData)));
        Percentile_95_AnalysedData{i}{j} = TempData(round(0.95*length(TempData)));
        IQRAnalysedData{i}{j} = iqr(TempData);
        
        Threshold{i}{j} = median(TempData) + OutlierThreshold*iqr(TempData);
        
        
        for k = 1:length(AnalysedData{i}{j}.UnDirBout.MaxBoutSeqMatch),
            [Pks{i}{j}{k}, Locs{i}{j}{k}] = findpeaks(AnalysedData{i}{j}.UnDirBout.MaxBoutSeqMatch{k}, 'MINPEAKHEIGHT', Threshold{i}{j});
        end
        for k = 1:length(AnalysedData{i}{j}.DirBout.MaxBoutSeqMatch),
          [Pks{i}{j}{k}, Locs{i}{j}{k}] = findpeaks(AnalysedData{i}{j}.DirBout.MaxBoutSeqMatch{k}, 'MINPEAKHEIGHT', Threshold{i}{j});
        end
    end
end

disp('Finished');

