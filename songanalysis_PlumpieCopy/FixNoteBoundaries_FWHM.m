function [TimeNow, LogFileName] = FixNoteBoundaries_FWHM(DataDir, RecFileDir, NoteDir, NoteFileList, FileType, PlotOption)

PresentDir = pwd;
!mkdir FixedNoteFiles_FWHM

PreDur = 0.01; % seconds
PostDur = 0.01; % seconds
LargeShiftAmount = 0.025; % seconds

Fid = fopen(NoteFileList, 'r');
TempNoteFiles = textscan(Fid, '%s', 'DeLimiter', '\n');
TempNoteFiles = TempNoteFiles{1};
fclose(Fid);

for i = 1:length(TempNoteFiles),
    NoteFiles(i).name = [TempNoteFiles{i}, '.not.mat'];
end

TimeNow = datestr(now, 'mmddyyHHMMSS');
LogFileName = [NoteFileList,'.FixNoteBoundaries_FWHMLog.', TimeNow, '.txt'];
LogFid = fopen(LogFileName, 'w');

fprintf(LogFid, 'Log file for the shifts of note boundaries fixed to ensure that all the note onsets and offsets are done using a threshold set to half the maximum of the syllable segmented with uisonganal_rr_fisher_threshold\n');
fprintf(LogFid, 'Positive shifts indicate that the new boundary was shifted forward and negative shifts indicate the new boundary was shifted back in time\n');

fprintf(LogFid, '\nData directory is %s\n', DataDir);
fprintf(LogFid, 'Rec file directory is %s\n', RecFileDir);
fprintf(LogFid, 'Note file directory is %s\n', NoteDir);
fprintf(LogFid, 'Song file list is %s\n', NoteFileList);
fprintf(LogFid, 'File type is %s\n', FileType);
fprintf(LogFid, 'Plot option is %s\n', PlotOption);
fprintf(LogFid, '\nPre duration to be considered is %g seconds\n', PreDur);
fprintf(LogFid, 'Post duration to be considered is %g seconds\n\n', PostDur);

SyllNo = 0;
LargeShiftSylls = 0;
NoShiftSylls = 0;

for i = 1:length(NoteFiles),
    Index = strfind(NoteFiles(i).name, '.not.mat');
    SongFile = NoteFiles(i).name(1:Index-1);
    load([NoteDir, '/', NoteFiles(i).name]);

    fprintf(LogFid, 'Song file name: %s\n', SongFile);
    fprintf(LogFid, '======================================================================\n', SongFile);
    fprintf(LogFid, 'Original threshold used by uisonganal_rr_fisher_threshold is %g\n\n', threshold);
    
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

    Song = resample(Song, 44100, Fs);
    Fs = 44100;
    Time = (1:1:length(Song))/Fs;
    
    FFTWinSize = 0.008;
    WinSize = round(FFTWinSize * Fs);
    WinOverlap = round(WinSize * 0.9);
    [S, F, T, P] = spectrogram(Song, hamming(WinSize), WinOverlap, WinSize, Fs);
    Freq = find((F >= 300) & (F <= 8000));
    smooth = 10*log10(sum((S(Freq,:)).*conj((S(Freq,:)))));
    smooth = spline(T, smooth, Time);
    
    Time = Time*1000;
       
    SmoothBaseline = sort(smooth);
    SmoothBaseline = SmoothBaseline(round(0.15*length(SmoothBaseline)));
    
    for j = 1:length(onsets),
        if (labels(j) == '0')
            continue;
        end
        fprintf(LogFid, 'Syllable %c:', labels(j));
        
        if ((onsets(j) - round(PreDur*Fs)) < Time(1))
            fprintf(LogFid, ' not enough pre time, so not shifted\n');
            NoShiftSylls = NoShiftSylls + 1;
            continue;
        end
        
        if ((offsets(j) + round(PostDur*Fs)) > Time(end))
            fprintf(LogFid, ' not enough post time, so not shifted\n');
            NoShiftSylls = NoShiftSylls + 1;
            continue;
        end
        
        StartIndex = find(Time >= onsets(j), 1, 'first');
        EndIndex = find(Time >= offsets(j), 1, 'first');
        
        MaxVal = max(smooth(StartIndex:EndIndex));
        
        Threshold = (MaxVal + SmoothBaseline)/2;
        
        fprintf(LogFid, ' threshold is %g:', Threshold);
        
        NewStartIndex = find(smooth((StartIndex - round(PreDur*Fs)):(StartIndex + round(PostDur*Fs))) >= Threshold, 1, 'first');
        NewStartIndex = NewStartIndex + (StartIndex - round(PreDur*Fs)) - 1;
        if (isempty(NewStartIndex))
            NewStartIndex = find(smooth(StartIndex:EndIndex) >= Threshold, 1, 'first');
            NewStartIndex = NewStartIndex + StartIndex - 1;
        end
                
        NewEndIndex = find(smooth((EndIndex - round(PreDur*Fs)):(EndIndex + round(PostDur*Fs))) <= Threshold, 1, 'first');
        NewEndIndex = NewEndIndex + (EndIndex - round(PreDur*Fs)) - 1;
        if (isempty(NewEndIndex))
            NewEndIndex = find(smooth(EndIndex:end) <= Threshold, 1, 'first');
            NewEndIndex = NewEndIndex + EndIndex - 1;
        end
        
        SyllNo = SyllNo + 1;
        Label(SyllNo) = labels(j);
        ShiftFlag = 0;
        OnsetShift(SyllNo) = Time(NewStartIndex) - onsets(j);
        if (abs(OnsetShift(SyllNo)) > (LargeShiftAmount*1000))
            NewStartIndex = StartIndex;
            OnsetShift(SyllNo) = -1000;
            fprintf(LogFid, ' onset shift too large, so did not shift: ');
            ShiftFlag = 1;
            LargeShiftSylls = LargeShiftSylls + 1;
        else
            fprintf(LogFid, ' onset shift %g ms: ', OnsetShift(SyllNo));
        end
                
        OffsetShift(SyllNo) = Time(NewEndIndex) - offsets(j);
        if (abs(OffsetShift(SyllNo)) > (LargeShiftAmount * 1000))
            NewEndIndex = EndIndex;
            OffsetShift(SyllNo) = -1000;
            fprintf(LogFid, ' offset shift too large, so did not shift: ');
            if (ShiftFlag == 0)
                LargeShiftSylls = LargeShiftSylls + 1;
            end
        else
            fprintf(LogFid, ' offset shift %g ms\n', OffsetShift(SyllNo));
        end
        
        onsets(j) = Time(NewStartIndex);
        offsets(j) = Time(NewEndIndex);
        
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
        for k = 1:length(onsets),
            plot([onsets(k)/1000 onsets(k)/1000 offsets(k)/1000 offsets(k)/1000], [0 7000 7000 0], 'k');
        end
        uiwait(gcf);
    end
    if (ispc)
        save(['FixedNoteFiles_FWHM\', NoteFiles(i).name], 'labels', 'Fs', 'offsets', 'onsets', 'TimeNow');    
    else
        save(['FixedNoteFiles_FWHM/', NoteFiles(i).name], 'labels', 'Fs', 'offsets', 'onsets', 'TimeNow');
    end
    fprintf(LogFid, '\n');
end
fprintf(LogFid, '======================================================================\n', SongFile);
fprintf(LogFid, '\nTotal # of syllables = %i\n', NoShiftSylls + SyllNo);
fprintf(LogFid, '# of syllables that were not shifted because there was not enough time pre or post = %i\n', NoShiftSylls);
fprintf(LogFid, '# of syllables that were not shifted because the shift was larger than %g ms = %i\n', LargeShiftAmount*1000, LargeShiftSylls);
fprintf(LogFid, '======================================================================\n', SongFile);

UniqueLabels = unique(Label);

fprintf(LogFid, '\nIndividual Syllable Shift Statistics\n\n');
for i = 1:length(UniqueLabels),
    Matches = find(Label == UniqueLabels(i));
    SyllOnsetShifts = OnsetShift(Matches);
    SyllOnsetShifts = SyllOnsetShifts(find(SyllOnsetShifts > -1000));
    SyllOnsetShifts = sort(SyllOnsetShifts);
    
    SyllOffsetShifts = OffsetShift(Matches);
    SyllOffsetShifts = SyllOffsetShifts(find(SyllOffsetShifts > -1000));
    SyllOffsetShifts = sort(SyllOffsetShifts);
    
    fprintf(LogFid,'Syllable %c:\n', UniqueLabels(i));
    fprintf(LogFid, 'Onset Shift - %i occurences\n', length(SyllOnsetShifts));
    fprintf(LogFid, 'Min = %g ms\n', min(SyllOnsetShifts));
    if (length(SyllOnsetShifts) >= 10)
        fprintf(LogFid, '5th percentile = %g ms\n', (SyllOnsetShifts(round(0.05*length(SyllOnsetShifts)))));
    end
    fprintf(LogFid, 'Mean = %g ms\n', mean(SyllOnsetShifts));
    fprintf(LogFid, 'Median = %g ms\n', median(SyllOnsetShifts));
    
    if (length(SyllOnsetShifts) >= 10)
        fprintf(LogFid, '95th percentile = %g ms\n', (SyllOnsetShifts(round(0.95*length(SyllOnsetShifts)))));
    end
    fprintf(LogFid, 'Max = %g ms\n\n', max(SyllOnsetShifts));
    
    fprintf(LogFid, 'Offset Shift - %i occurences\n',  length(SyllOffsetShifts));
    fprintf(LogFid, 'Min = %g ms\n', min(SyllOffsetShifts));
    if (length(SyllOffsetShifts) >= 10)
        fprintf(LogFid, '5th percentile = %g ms\n', (SyllOffsetShifts(round(0.05*length(SyllOffsetShifts)))));
    end
    fprintf(LogFid, 'Mean = %g ms\n', mean(SyllOffsetShifts));
    fprintf(LogFid, 'Median = %g ms\n', median(SyllOffsetShifts));
    if (length(SyllOffsetShifts) >= 10)
        fprintf(LogFid, '95th percentile = %g ms\n', (SyllOffsetShifts(round(0.95*length(SyllOffsetShifts)))));
    end
    fprintf(LogFid, 'Max = %g ms\n\n', max(SyllOffsetShifts));
    
%    figure;
%    subplot(2,1,1);

    Edges = min(SyllOnsetShifts)-2:0.25:max(SyllOnsetShifts)+2;
%    SyllOnsetShiftBar = bar(Edges, histc(SyllOnsetShifts,Edges)*100/sum(histc(SyllOnsetShifts,Edges)), 'histc');
%    set(SyllOnsetShiftBar, 'FaceColor', 'none', 'EdgeColor', 'k');
%    hold on;
%    axis tight;
%    TempAxis = axis;
%    plot([min(SyllOnsetShifts) min(SyllOnsetShifts)], [TempAxis(3) TempAxis(4)], 'r:', 'LineWidth', 2);
%    plot([max(SyllOnsetShifts) max(SyllOnsetShifts)], [TempAxis(3) TempAxis(4)], 'r:', 'LineWidth', 2);
%    plot([median(SyllOnsetShifts) median(SyllOnsetShifts)], [TempAxis(3) TempAxis(4)], 'b--', 'LineWidth', 2);
%    plot([mean(SyllOnsetShifts) mean(SyllOnsetShifts)], [TempAxis(3) TempAxis(4)], 'b--', 'LineWidth', 2);

    SyllOnsetShifts = sort(SyllOnsetShifts);
    if (length(SyllOnsetShifts) >= 10)
%        plot([SyllOnsetShifts(round(0.95*length(SyllOnsetShifts))) SyllOnsetShifts(round(0.95*length(SyllOnsetShifts)))], [TempAxis(3) TempAxis(4)], 'm--', 'LineWidth', 2);
%        plot([SyllOnsetShifts(round(0.05*length(SyllOnsetShifts))) SyllOnsetShifts(round(0.05*length(SyllOnsetShifts)))], [TempAxis(3) TempAxis(4)], 'm--', 'LineWidth', 2);
    end
    
%    xlabel('Time (ms)');
%    ylabel('%');
%    title(['Shifts in onsets for ', num2str(length(SyllOnsetShifts)), ' occurences of syllable ', UniqueLabels(i)]);

%    subplot(2,1,2);
%    Edges = min(SyllOffsetShifts)-2:1:max(SyllOffsetShifts)+2;
%    SyllOffsetShiftBar = bar(Edges, histc(SyllOffsetShifts,Edges)*100/sum(histc(SyllOffsetShifts,Edges)), 'histc');
%    set(SyllOffsetShiftBar, 'FaceColor', 'none', 'EdgeColor', 'k');
%    hold on;
%    axis tight;
%    TempAxis = axis;
%    plot([min(SyllOffsetShifts) min(SyllOffsetShifts)], [TempAxis(3) TempAxis(4)], 'r:', 'LineWidth', 2);
%    plot([max(SyllOffsetShifts) max(SyllOffsetShifts)], [TempAxis(3) TempAxis(4)], 'r:', 'LineWidth', 2);
%    plot([median(SyllOffsetShifts) median(SyllOffsetShifts)], [TempAxis(3) TempAxis(4)], 'b--', 'LineWidth', 2);
%    plot([mean(SyllOffsetShifts) mean(SyllOffsetShifts)], [TempAxis(3) TempAxis(4)], 'b--', 'LineWidth', 2);

    SyllOffsetShifts = sort(SyllOffsetShifts);
    if (length(SyllOffsetShifts) >= 10)
%        plot([SyllOffsetShifts(round(0.95*length(SyllOffsetShifts))) SyllOffsetShifts(round(0.95*length(SyllOffsetShifts)))], [TempAxis(3) TempAxis(4)], 'm--', 'LineWidth', 2);
%        plot([SyllOffsetShifts(round(0.05*length(SyllOffsetShifts))) SyllOffsetShifts(round(0.05*length(SyllOffsetShifts)))], [TempAxis(3) TempAxis(4)], 'm--', 'LineWidth', 2);
    end

 %   xlabel('Time (ms)');
 %   ylabel('%');
 %   title(['Shifts in offsets for ', num2str(length(SyllOffsetShifts)), ' occurences of syllable ', UniqueLabels(i)]);
 %   saveas(gcf, [NoteFileList, '.FixNoteBoundaries_FWHM_shifts_summary_', TimeNow, '.Syllable', UniqueLabels(i), '.fig']);
 %   saveas(gcf, [NoteFileList, '.FixNoteBoundaries_FWHM_shifts_summary_', TimeNow, '.Syllable', UniqueLabels(i), '.png'], 'png');
    c;
end

%figure;
%subplot(2,1,1);
TempOnsetShift = OnsetShift;
OnsetShift = OnsetShift(find(OnsetShift > -1000));

Edges = min(OnsetShift)-2:0.25:max(OnsetShift)+2;
%OnsetShiftBar = bar(Edges, histc(OnsetShift,Edges)*100/sum(histc(OnsetShift,Edges)), 'histc');
%set(OnsetShiftBar, 'FaceColor', 'none', 'EdgeColor', 'k');
%hold on;
%axis tight;
%TempAxis = axis;
%plot([min(OnsetShift) min(OnsetShift)], [TempAxis(3) TempAxis(4)], 'r:', 'LineWidth', 2);
%plot([max(OnsetShift) max(OnsetShift)], [TempAxis(3) TempAxis(4)], 'r:', 'LineWidth', 2);
%plot([median(OnsetShift) median(OnsetShift)], [TempAxis(3) TempAxis(4)], 'b--', 'LineWidth', 2);
%plot([mean(OnsetShift) mean(OnsetShift)], [TempAxis(3) TempAxis(4)], 'b--', 'LineWidth', 2);

SortedOnsetShifts = sort(OnsetShift);
%plot([SortedOnsetShifts(round(0.95*length(SortedOnsetShifts))) SortedOnsetShifts(round(0.95*length(SortedOnsetShifts)))], [TempAxis(3) TempAxis(4)], 'm--', 'LineWidth', 2);
%plot([SortedOnsetShifts(round(0.05*length(SortedOnsetShifts))) SortedOnsetShifts(round(0.05*length(SortedOnsetShifts)))], [TempAxis(3) TempAxis(4)], 'm--', 'LineWidth', 2);
%xlabel('Time (ms)');
%ylabel('%');
%title(['Shifts in onsets for ', num2str(length(OnsetShift)), ' syllables']);

fprintf(LogFid, '======================================================================\n', SongFile);
fprintf(LogFid, 'Overall shift statistics for all syllables\n\n');

fprintf(LogFid, 'Onset Shift - %i occurences\n', length(OnsetShift));
fprintf(LogFid, 'Min = %g ms\n', min(OnsetShift));
fprintf(LogFid, '5th percentile = %g ms\n', (SortedOnsetShifts(round(0.05*length(SortedOnsetShifts)))));
fprintf(LogFid, 'Mean = %g ms\n', mean(OnsetShift));
fprintf(LogFid, 'Median = %g ms\n', median(OnsetShift));
fprintf(LogFid, '95th percentile = %g ms\n', (SortedOnsetShifts(round(0.95*length(SortedOnsetShifts)))));
fprintf(LogFid, 'Max = %g ms\n\n', max(OnsetShift));

%subplot(2,1,2);
%TempOffsetShift = OffsetShift;
%OffsetShift = OffsetShift(find(OffsetShift > -1000));

%Edges = min(OffsetShift)-2:1:max(OffsetShift)+2;
%OffsetShiftBar = bar(Edges, histc(OffsetShift,Edges)*100/sum(histc(OffsetShift,Edges)), 'histc');
%set(OffsetShiftBar, 'FaceColor', 'none', 'EdgeColor', 'k');
%hold on;
%axis tight;
%TempAxis = axis;
%plot([min(OffsetShift) min(OffsetShift)], [TempAxis(3) TempAxis(4)], 'r:', 'LineWidth', 2);
%plot([max(OffsetShift) max(OffsetShift)], [TempAxis(3) TempAxis(4)], 'r:', 'LineWidth', 2);
%plot([median(OffsetShift) median(OffsetShift)], [TempAxis(3) TempAxis(4)], 'b--', 'LineWidth', 2);
%plot([mean(OffsetShift) mean(OffsetShift)], [TempAxis(3) TempAxis(4)], 'b--', 'LineWidth', 2);

SortedOffsetShifts = sort(OffsetShift);
%plot([SortedOffsetShifts(round(0.95*length(SortedOffsetShifts))) SortedOffsetShifts(round(0.95*length(SortedOffsetShifts)))], [TempAxis(3) TempAxis(4)], 'm--', 'LineWidth', 2);
%plot([SortedOffsetShifts(round(0.05*length(SortedOffsetShifts))) SortedOffsetShifts(round(0.05*length(SortedOffsetShifts)))], [TempAxis(3) TempAxis(4)], 'm--', 'LineWidth', 2);

xlabel('Time (ms)');
ylabel('%');
%title(['Shifts in offsets for ', num2str(length(OffsetShift)), ' syllables']);
%saveas(gcf, [NoteFileList, '.FixNoteBoundaries_FWHM_shifts_summary_', TimeNow, '.Overall.fig']);
%saveas(gcf, [NoteFileList, '.FixNoteBoundaries_FWHM_shifts_summary_', TimeNow, '.Overall.png'], 'png');
c;
    
fprintf(LogFid, 'Offset Shift - %i occurences\n', length(OffsetShift));
fprintf(LogFid, 'Min = %g ms\n', min(OffsetShift));
fprintf(LogFid, '5th percentile = %g ms\n', (SortedOffsetShifts(round(0.05*length(SortedOffsetShifts)))));
fprintf(LogFid, 'Mean = %g ms\n', mean(OffsetShift));
fprintf(LogFid, 'Median = %g ms\n', median(OffsetShift));
fprintf(LogFid, '95th percentile = %g ms\n', (SortedOffsetShifts(round(0.95*length(SortedOffsetShifts)))));
fprintf(LogFid, 'Max = %g ms\n\n', max(OffsetShift));

fclose(LogFid);
disp('Finished saving fixed note files');
