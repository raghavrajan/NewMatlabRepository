function [] = AnalyseSyllableSequence(Dir, BirdName, Motif, FilesToAnalyse)

PreWindow = 0.5; % in sec
PostWindow = 0.5; % in sec

cd(Dir);

for i = 1:length(Motif),
    Files = dir(['*Syll', Motif(i), '*.mat']);
%     RandomFileName = Files(end).name;
%     disp(RandomFileName);
%     if (iscell(BirdName))
%         [MeanBoutPeaks, StdBoutPeaks, RandPeakVals, RandTotalTime] = CalculateMeanRandomMatches(RandomFileName, 'k', BirdName, 'off', ThresholdMultiplier(i));
%     else
%         [MeanBoutPeaks, StdBoutPeaks, RandPeakVals, RandTotalTime] = CalculateMeanRandomMatches(RandomFileName, 'k', {BirdName}, 'off', ThresholdMultiplier(i));
%     end
%     Threshold(i) = MeanBoutPeaks + ThresholdMultiplier(i)*StdBoutPeaks;
    
    for j = 1:length(Files) - 1,
        Sylls{i}{j} = load(Files(j).name);
    end
end

for i = 1:length(Motif),
    disp(['Undir syllable is ', Motif(i)]);
    figure;
    TempPeaks = [];
    for TempCounter = 1:length(Sylls{i}{1}.UnDirBout.Peaks),
        TempPeaks = [TempPeaks Sylls{i}{1}.UnDirBout.Peaks{TempCounter}];
    end
    Edges = linspace(0,1.1*max(TempPeaks));
    plot(Edges, histc(TempPeaks, Edges)/sum(histc(TempPeaks, Edges)), 'b');
    hold on;
    
    TempPeaks = [];
    for TempCounter = 1:length(Sylls{i}{1}.DirBout.Peaks),
        TempPeaks = [TempPeaks Sylls{i}{1}.DirBout.Peaks{TempCounter}];
    end
    Edges = linspace(0,1.1*max(TempPeaks));
    plot(Edges, histc(TempPeaks, Edges)/sum(histc(TempPeaks, Edges)), 'r');
    
    uiwait(gcf);
    Threshold(i) = input('What is the threshold you want to use for this syllable: ');
end
Colours = ['rgbkcm'];
for i = 1:length(Motif),
    disp(['Undir Syllable is ', Motif(i), ' and the threshold is ', num2str(Threshold(i))]);
    for FileNo = FilesToAnalyse,
        disp(['Day # ', num2str(FileNo)]);
        RemovalIndices = 0;
        NoofMatches = 0;
        for k = 1:length(Motif),
            UnDirMatches{i}{FileNo}{k} = [];
        end
        for j = 1:length(Sylls{i}{FileNo}.UnDirBout.Peaks);
            Locs = find(Sylls{i}{FileNo}.UnDirBout.Peaks{j} > Threshold(i));
            if (size(Locs,1) < size(Locs,2))
                Locs = Locs';
            end
            LocTimes = Sylls{i}{FileNo}.UnDirBout.PeakTimes{j}(Locs);
            
            if (~isempty(Locs))
                RemovalIndices = RemovalIndices + length(find((LocTimes < PreWindow) & (LocTimes > (Sylls{i}{FileNo}.UnDirBout.FileLength{j} - PostWindow))));
                TempLocTimes = LocTimes((LocTimes > PreWindow) & (LocTimes < (Sylls{i}{FileNo}.UnDirBout.FileLength{j} - PostWindow)));
                if (~isempty(TempLocTimes))
                    for LocNo = 1:length(TempLocTimes)
                        NoofMatches = NoofMatches + 1;
                        for k = 1:length(Motif),
                            SyllableMatches = Sylls{k}{FileNo}.UnDirBout.Peaks{j}(find((Sylls{k}{FileNo}.UnDirBout.PeakTimes{j} > (TempLocTimes(LocNo) - PreWindow)) & (Sylls{k}{FileNo}.UnDirBout.PeakTimes{j} < (TempLocTimes(LocNo) + PostWindow))));
                            SyllableMatchTimes = Sylls{k}{FileNo}.UnDirBout.PeakTimes{j}(find((Sylls{k}{FileNo}.UnDirBout.PeakTimes{j} > (TempLocTimes(LocNo) - PreWindow)) & (Sylls{k}{FileNo}.UnDirBout.PeakTimes{j} < (TempLocTimes(LocNo) + PostWindow))));
                            ActualSyllableMatches = (find(SyllableMatches > Threshold(k)));
                            if (size(SyllableMatches,1) < size(SyllableMatches,2))
                                SyllableMatches = SyllableMatches';
                                SyllableMatchTimes = SyllableMatchTimes';
                            end
                            if (~isempty(ActualSyllableMatches))
                                UnDirMatches{i}{FileNo}{k} = [UnDirMatches{i}{FileNo}{k}; [(SyllableMatchTimes(ActualSyllableMatches) - TempLocTimes(LocNo)) SyllableMatches(ActualSyllableMatches)]];
                            end
                        end
                    end
                end
            end
        end
        disp(['Used ', num2str(NoofMatches), ' syllables']);
        disp(['Removed ', num2str(RemovalIndices), ' syllables since they did not have enough data before or after']);
    end
end

for i = 1:length(Motif),
    disp(['Undir Syllable is ', Motif(i), ' and the threshold is ', num2str(Threshold(i))]);
    for FileNo = FilesToAnalyse,
        disp(['Day # ', num2str(FileNo)]);
        RemovalIndices = 0;
        NoofMatches = 0;
        for k = 1:length(Motif),
            DirMatches{i}{FileNo}{k} = [];
        end
        for j = 1:length(Sylls{i}{FileNo}.DirBout.Peaks);
            Locs = find(Sylls{i}{FileNo}.DirBout.Peaks{j} > Threshold(i));
            if (size(Locs,1) < size(Locs,2))
                Locs = Locs';
            end
            LocTimes = Sylls{i}{FileNo}.DirBout.PeakTimes{j}(Locs);
            
            if (~isempty(Locs))
                RemovalIndices = RemovalIndices + length(find((LocTimes < PreWindow) & (LocTimes > (Sylls{i}{FileNo}.DirBout.FileLength{j} - PostWindow))));
                TempLocTimes = LocTimes((LocTimes > PreWindow) & (LocTimes < (Sylls{i}{FileNo}.DirBout.FileLength{j} - PostWindow)));
                if (~isempty(TempLocTimes))
                    for LocNo = 1:length(TempLocTimes)
                        NoofMatches = NoofMatches + 1;
                        for k = 1:length(Motif),
                            SyllableMatches = Sylls{k}{FileNo}.DirBout.Peaks{j}(find((Sylls{k}{FileNo}.DirBout.PeakTimes{j} > (TempLocTimes(LocNo) - PreWindow)) & (Sylls{k}{FileNo}.DirBout.PeakTimes{j} < (TempLocTimes(LocNo) + PostWindow))));
                            SyllableMatchTimes = Sylls{k}{FileNo}.DirBout.PeakTimes{j}(find((Sylls{k}{FileNo}.DirBout.PeakTimes{j} > (TempLocTimes(LocNo) - PreWindow)) & (Sylls{k}{FileNo}.DirBout.PeakTimes{j} < (TempLocTimes(LocNo) + PostWindow))));
                            ActualSyllableMatches = (find(SyllableMatches > Threshold(k)));
                            if (size(SyllableMatches,1) < size(SyllableMatches,2))
                                SyllableMatches = SyllableMatches';
                                SyllableMatchTimes = SyllableMatchTimes';
                            end
                            if (~isempty(ActualSyllableMatches))
                                DirMatches{i}{FileNo}{k} = [DirMatches{i}{FileNo}{k}; [(SyllableMatchTimes(ActualSyllableMatches) - TempLocTimes(LocNo)) SyllableMatches(ActualSyllableMatches)]];
                            end
                        end
                    end
                end
            end
        end
        disp(['Used ', num2str(NoofMatches), ' syllables']);
        disp(['Removed ', num2str(RemovalIndices), ' syllables since they did not have enough data before or after']);
    end
end

Figure = figure;
set(gcf, 'Color', 'w');
for i = 1:length(Motif),
    
    if (i == 1)
        NextSyll = i+1;
        PrevSyll = length(Motif);
    else
        if (i == length(Motif))
            NextSyll = 1;
            PrevSyll = i - 1;
        else
            NextSyll = i+1;
            PrevSyll = i-1;
        end
    end
    
    subplot(length(Motif),2,(i-1)*2 + 1);
    hold on;
    for FileNo = FilesToAnalyse,
        UnDirBar{i}{FileNo}{1} = bar((FileNo*2-1),(length(find(UnDirMatches{i}{FileNo}{PrevSyll}(:,1) < 0))/size(UnDirMatches{i}{FileNo}{i},1)));
        set(UnDirBar{i}{FileNo}{1}, 'FaceColor', 'none', 'EdgeColor', 'b');
        DirBar{i}{FileNo}{1} = bar((FileNo*2),(length(find(DirMatches{i}{FileNo}{PrevSyll}(:,1) < 0))/size(DirMatches{i}{FileNo}{i},1)));
        set(DirBar{i}{FileNo}{1}, 'FaceColor', 'none', 'EdgeColor', 'r');
    end
    title([Motif(i), ' preceded by ', Motif(PrevSyll)], 'FontSize', 12, 'FontName', 'Arial');
    ylabel('Fraction of total syllables', 'FontSize', 8, 'FontName', 'Arial');
    
    axis([0.4 (FileNo*2 + 0.6) 0 1.1]);
    
    subplot(length(Motif),2,(i-1)*2 + 2);
    hold on;
    for FileNo = FilesToAnalyse,
        UnDirBar{i}{FileNo}{2} = bar((FileNo*2-1),(length(find(UnDirMatches{i}{FileNo}{NextSyll}(:,1) > 0))/size(UnDirMatches{i}{FileNo}{i},1)));
        set(UnDirBar{i}{FileNo}{2}, 'FaceColor', 'none', 'EdgeColor', 'b');
        DirBar{i}{FileNo}{2} = bar((FileNo*2),(length(find(DirMatches{i}{FileNo}{NextSyll}(:,1) > 0))/size(DirMatches{i}{FileNo}{i},1)));
        set(DirBar{i}{FileNo}{2}, 'FaceColor', 'none', 'EdgeColor', 'r');
    end
    title([Motif(i), ' followed by ', Motif(NextSyll)], 'FontSize', 12, 'FontName', 'Arial');
    ylabel('Fraction of total syllables', 'FontSize', 8, 'FontName', 'Arial');
    axis([0.4 (FileNo*2 + 0.6) 0 1.1]);
end