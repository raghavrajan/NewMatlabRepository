function [MatchLabels] = AutoSegmentandLabelSongs(DataDir, RecFileDir, SongFileList, FileType, Templates, Labels, StretchValues, PlotOption)

min_int = 0.005; % in sec
min_dur = 0.01; % in sec

PresentDir = pwd;
!mkdir AutoLabeledNoteFiles

Fid = fopen(SongFileList, 'r');
TempFiles = textscan(Fid, '%s', 'DeLimiter', '\n');
SongFiles = TempFiles{1};
fclose(Fid);

TimeNow = datestr(now, 'mmddyyHHMMSS');
%LogFileName = [SongFileList,'.AutoSegmentandLabelSongs.log.', TimeNow, '.txt'];
%LogFid = fopen(LogFileName, 'w');

%fprintf(LogFid, 'Log file for the shifts of note boundaries fixed to ensure that all the note onsets and offsets are done using a threshold set to half the maximum of the syllable segmented with uisonganal_rr_fisher_threshold\n');
%fprintf(LogFid, 'Positive shifts indicate that the new boundary was shifted forward and negative shifts indicate the new boundary was shifted back in time\n');

% fprintf(LogFid, '\nData directory is %s\n', DataDir);
% fprintf(LogFid, 'Rec file directory is %s\n', RecFileDir);
% fprintf(LogFid, 'Note file directory is %s\n', NoteDir);
% fprintf(LogFid, 'Song file list is %s\n', NoteFileList);
% fprintf(LogFid, 'File type is %s\n', FileType);
% fprintf(LogFid, 'Plot option is %s\n', PlotOption);
% fprintf(LogFid, '\nPre duration to be considered is %g seconds\n', PreDur);
% fprintf(LogFid, 'Post duration to be considered is %g seconds\n\n', PostDur);

for i = 1:length(SongFiles),
    SongFile = SongFiles{i};
    SlashIndex = find((SongFile == '/') | (SongFile == '\'));
    if (~isempty(SlashIndex))
        SongFile = SongFile(SlashIndex(end)+1:end);
    end
    if (strfind(FileType, 'okrank'))
        [Song, Fs] = SSAReadOKrankData(DataDir, RecFileDir, SongFile, 1);
    else
        disp(SongFile);
        if (strfind(FileType, 'wav'))
            cd(DataDir);
            [Song, Fs] = wavread(SongFile);
            cd(PresentDir);
        else
            if (strfind(FileType, 'obs'));
                [Song, Fs] = SSASoundIn(DataDir, RecFileDir, SongFile, 'obs0r');
                Song = Song * 5/32768;
            end
        end
    end

    Time = (1:1:length(Song))/Fs;
    
    FFTWinSize = 0.008;
    WinSize = round(FFTWinSize * Fs);
    WinOverlap = round(WinSize * 0.9);
    [S, F, T, P] = spectrogram(Song, hamming(WinSize), WinOverlap, WinSize, Fs);
    Freq = find((F >= 300) & (F <= 8000));
    smooth = 10*log10(sum((S(Freq,:)).*conj((S(Freq,:)))));
    smooth = spline(T, smooth, Time);
    
    Obj = gmdistribution.fit(smooth', 2);
    
    threshold = mean(Obj.mu);
    
    [onsets, offsets] = segment(smooth, Fs, min_int*1000, min_dur*1000, threshold);
    onsets = onsets/1000;
    offsets = offsets/1000;
    
    clear SyllMatch;
    clear labels;
    for j = 1:length(onsets),
    
        fprintf('%i', j);
        FFTWinOverlap = 0.5;
        Syllable = Song(find((Time >= (onsets(j)-0.020)) & (Time <= (offsets(j) + 0.020))));
        WinSize = round(FFTWinSize * Fs);
        WinOverlap = round(FFTWinOverlap * WinSize);
        [S1, F, T, P] = spectrogram(Syllable, hamming(WinSize), WinOverlap, WinSize, Fs);
        Freq1 = find((F >= 860) & (F <= 8600));
        Power = log10(sum(S1(Freq1,:).*conj(S1(Freq1,:))));
        S = log10(abs(S1(Freq1,:)));
        TempS = S;

%         figure;
%         contourf(TempS);
%         uiwait(gcf);

        for TemplateNo = 1:length(Templates),
            MotifTemplate = Templates{TemplateNo}.MotifTemplate;
            clear TempMatch;
            TemplateTime = 1:1:size(MotifTemplate, 2);
            StretchIndex = 0;
            for Stretch = [StretchValues],
                fprintf('>');
                StretchIndex = StretchIndex + 1;
                STemplateTime = linspace(TemplateTime(1), TemplateTime(end), (size(MotifTemplate,2) * (1 + Stretch/100)));
                clear WMotif;
                WMotif = zeros(size(MotifTemplate,1), length(STemplateTime));
                for WMotifRow = 1:size(MotifTemplate, 1);
                    WMotif(WMotifRow,:) = interp1(TemplateTime, MotifTemplate(WMotifRow,:), STemplateTime);
                end
        
                if (size(TempS,2) < size(WMotif,2))
                    Match = [];
                else
                    WinMean = zeros((size(TempS,2) - size(WMotif,2) + 1), 1);
                    WinSTD = zeros((size(TempS,2) - size(WMotif,2) + 1), 1);

                    TempMeanSTD = CalculateMeanSTDforSpectralMatch(TempS(1:size(TempS,1)*size(TempS,2)), size(WMotif,1)*size(WMotif,2), (size(TempS,2) - size(WMotif,2) + 1), size(WMotif,1));

                    WinMean = TempMeanSTD(1:length(TempMeanSTD)/2);
                    WinSTD = TempMeanSTD((length(TempMeanSTD)/2 + 1):end);
                    [Match] = CalTemplateMatch(WMotif, TempS, WinMean, WinSTD);
                    Match = Match*size(WMotif,1)*size(WMotif,2);
                end
                if (~isempty(Match))
                    TempMatch(StretchIndex) = max(Match);
                else
                    TempMatch(StretchIndex) = NaN;
                end
            end
            SyllMatch{j}(TemplateNo) = max(TempMatch);
            fprintf('\t');
        end
        fprintf('\n');
        [MaxVal, MaxIndex] = max(SyllMatch{j});
        labels(j) = Labels(MaxIndex);
    end
    
    if (strfind(PlotOption, 'on'))
        if (strfind(FileType, 'wav'))
            cd (DataDir);
            if (ispc)
                PlotSpectrogram([DataDir, '\'], SongFile, FileType);
            else
                PlotSpectrogram([DataDir, '/'], SongFile, FileType);
            end
            cd (PresentDir);
        else
            if (ispc)
                PlotSpectrogram([DataDir, '\'], SongFile, FileType);
            else
                PlotSpectrogram([DataDir, '/'], SongFile, FileType);
            end
        end
        hold on;
        for j = 1:length(onsets),
            plot([onsets(j) onsets(j) offsets(j) offsets(j)], [0 7000 7000 0], 'k');
            text(mean([onsets(j) offsets(j)]), 7500, labels(j));
        end
        uiwait(gcf);
    end
    MatchLabels{i}.SongFile = SongFile;
    MatchLabels{i}.onsets = onsets;
    MatchLabels{i}.offsets = offsets;
    MatchLabels{i}.labels = labels;
    MatchLabels{i}.SyllMatch = SyllMatch;
end

disp('Finished labelling song files');
