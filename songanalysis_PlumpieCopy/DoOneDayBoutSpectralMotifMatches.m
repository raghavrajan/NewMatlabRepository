function [DirBout, UnDirBout] = DoOneDayBoutSpectralMotifMatches(PresentDirectory, MotifTemplate, SyllLabel, DirectoryName, DirFileList, UnDirFileList, FileType, DirBoutLimits, UnDirBoutLimits, StretchVals, GetRidNoise, PlotOption, MotifTemplateRawData, MotifTemplateRawDataFs)

cd(PresentDirectory);

Width = 0.005;
GaussianLen = 4;
Fs = 250;

XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * Fs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs) * (Width * Fs)));

for i = 1:length(MotifTemplate),
    cd(PresentDirectory);
%    UnDirBout{i} = BoutMotifSpectralMatch(DirectoryName, UnDirFileList, FileType, PlotOption, MotifTemplate{i}, StretchVals, GetRidNoise, UnDirBoutLimits, MotifTemplateRawData{i}, MotifTemplateRawDataFs{i});
    UnDirBout{i} = BoutMotifSpectralMatch(DirectoryName, UnDirFileList, FileType, PlotOption, MotifTemplate{i}, StretchVals, GetRidNoise, UnDirBoutLimits);
    cd(PresentDirectory);
%    DirBout{i} = BoutMotifSpectralMatch(DirectoryName, DirFileList, FileType, PlotOption, MotifTemplate{i}, StretchVals, GetRidNoise, DirBoutLimits, MotifTemplateRawData{i}, MotifTemplateRawDataFs{i});
    DirBout{i} = BoutMotifSpectralMatch(DirectoryName, DirFileList, FileType, PlotOption, MotifTemplate{i}, StretchVals, GetRidNoise, DirBoutLimits);
    cd(PresentDirectory);
end

for j = 1:length(UnDirBout{1}.MaxBoutSeqMatch),
    for i = 1:length(MotifTemplate),
        Durations(i,j) = length(UnDirBout{i}.MaxBoutSeqMatch{j});
    end
end

if (length(MotifTemplate) > 1)
    Durations = min(Durations);
end
for j = 1:length(UnDirBout{1}.MaxBoutSeqMatch),
    Match = zeros(length(MotifTemplate), Durations(j));
    for i = 1:length(MotifTemplate),
        Match(i,:) = UnDirBout{i}.MaxBoutSeqMatch{j}(1:Durations(j));
    end
    if (length(MotifTemplate) > 1)
        Match = max(Match);
    end
    %Match = conv(Match, GaussWin, 'same');
    %[Peaks, Locs] = findpeaks(Match, 'MINPEAKHEIGHT', 0);
    %RealUnDirBout.Peaks{j} = Peaks;
    %RealUnDirBout.PeakTimes{j} = UnDirBout{1}.T{j}(Locs);
    RealUnDirBout.MaxBoutSeqMatch{j} = Match;
    RealUnDirBout.T{j} = UnDirBout{1}.T{j}(1:length(Match));
    RealUnDirBout.FileName{j} = UnDirBout{1}.FileName{j};
    RealUnDirBout.MaxBoutSeqMatchVal{j} = UnDirBout{1}.MaxBoutSeqMatchVal{j};
    RealUnDirBout.FileLength{j} = UnDirBout{1}.FileLength{j};
end

clear Durations;
for j = 1:length(DirBout{1}.MaxBoutSeqMatch),
    for i = 1:length(MotifTemplate),
        Durations(i,j) = length(DirBout{i}.MaxBoutSeqMatch{j});
    end
end
if (length(MotifTemplate) > 1)
    Durations = min(Durations);
end

for j = 1:length(DirBout{1}.MaxBoutSeqMatch),
    Match = zeros(length(MotifTemplate), Durations(j));
    for i = 1:length(MotifTemplate),
        Match(i,:) = DirBout{i}.MaxBoutSeqMatch{j}(1:Durations(j));
    end
    if (length(MotifTemplate) > 1)
        Match = max(Match);
    end
    % Match = conv(Match, GaussWin, 'same');
%    [Peaks, Locs] = findpeaks(Match, 'MINPEAKHEIGHT', 0);
%    RealDirBout.Peaks{j} = Peaks;
%    RealDirBout.PeakTimes{j} = DirBout{1}.T{j}(Locs);
    RealDirBout.MaxBoutSeqMatch{j} = Match;
    RealDirBout.T{j} = DirBout{1}.T{j}(1:length(Match));
    RealDirBout.FileName{j} = DirBout{1}.FileName{j};
    RealDirBout.MaxBoutSeqMatchVal{j} = DirBout{1}.MaxBoutSeqMatchVal{j};
    RealDirBout.FileLength{j} = DirBout{1}.FileLength{j};
end

DirBout = RealDirBout;
UnDirBout = RealUnDirBout;

for i = 1:min([length(DirFileList) length(UnDirFileList)]);
    MatchFlag = strncmp(DirFileList, UnDirFileList, i);
    if (MatchFlag == 0)
        break;
    end
end

save([DirFileList(1:(i-2)), '.motifseqmatches', SyllLabel, '.mat'], 'UnDirBout', 'DirBout', 'DirectoryName', 'UnDirFileList', 'DirFileList', 'StretchVals');
% %================================================================%
