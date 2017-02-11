function [] = FixNoteBoundaries(DataDir, RecFileDir, NoteDir, NoteFileList, FileType, PlotOption)

PresentDir = pwd;
!mkdir FixedNoteFiles

ONMultiplier = 1/10;
OFFMultiplier = 1/12;
PreDur = 0.005; % seconds
PostDur = 0.005; % seconds

Fid = fopen(NoteFileList, 'r');
TempNoteFiles = textscan(Fid, '%s', 'DeLimiter', '\n');
TempNoteFiles = TempNoteFiles{1};
fclose(Fid);

for i = 1:length(TempNoteFiles),
    NoteFiles(i).name = [TempNoteFiles{i}, '.not.mat'];
end

for i = 1:length(NoteFiles),
    Index = strfind(NoteFiles(i).name, '.not.mat');
    SongFile = NoteFiles(i).name(1:Index-1);
    
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
    Time = (1:1:length(Song)-1)/Fs;
    
    FFTWinSize = 0.008;
    WinSize = round(FFTWinSize * Fs);
    WinOverlap = round(WinSize * 0.9);
    [S, F, T, P] = spectrogram(Song, hamming(WinSize), WinOverlap, WinSize, Fs);
    Freq = find((F >= 300) & (F <= 8000));
    smooth = 10*log10(sum((S(Freq,:)).*conj((S(Freq,:)))));
    Spect = 10*log10(S(Freq,:).*conj(S(Freq,:)));
    Deriv = diff(Spect, 1, 2);
    Deriv = sum(Deriv);
    %Deriv = diff(smooth);
    Window = ones(20,1);
    Window = Window/sum(Window);
    Deriv = conv(Deriv, Window, 'same');
    Deriv = conv(Deriv, Window, 'same');
    Deriv = spline(T(1:end-1),Deriv, Time);
    
    [pks, pklocs] = findpeaks(Deriv, 'MINPEAKHEIGHT', mean(Deriv) + std(Deriv));
    [troughs, troughlocs] = findpeaks(-Deriv, 'MINPEAKHEIGHT', -(mean(Deriv) - std(Deriv)));
        
    Time = Time*1000;
    
    Onsets = Time(pklocs);
    Offsets = Time(troughlocs);
    
    load([NoteDir, '/', NoteFiles(i).name]);
    
    Fs = 44100;
    
    Fid = fopen(['FixedNoteFiles/', SongFile, '.FixNoteBoundaries.log'], 'w');
 
    for j = 1:length(onsets),
    
        [MinVal, MinInd] = min(abs(Onsets - onsets(j)));
        
        if (abs(onsets(j) - Onsets(MinInd)) <= 15)
            %disp(['Shifted onset of Syllable ', labels(j), ' by ', num2str(onsets(j) - Onsets(MinInd)), ' ms']);
            fprintf(Fid, 'Shifted onset of Syllable %c by %g ms\n', labels(j), onsets(j) - Onsets(MinInd));
            onsets(j) = Onsets(MinInd);
        else
            disp(['Did not shift onset of Syllable ', labels(j)]);
            fprintf(Fid, 'Did not shift onset of Syllable %c\n', labels(j));
        end
        
        [MinVal, MinInd] = min(abs(Offsets - offsets(j)));
        if (abs(offsets(j) - Offsets(MinInd)) <= 15)    
            %disp(['Shifted offset of Syllable ', labels(j), ' by ', num2str(offsets(j) - Offsets(MinInd)), ' ms']);
            fprintf(Fid, 'Shifted offset of Syllable %c by %g ms\n', labels(j), offsets(j) - Offsets(MinInd));
            offsets(j) = Offsets(MinInd);
        else
            disp(['Did not shift offset of Syllable ', labels(j)]);
            fprintf(Fid, 'Did not shift offset of Syllable %c\n', labels(j));
        end
    end
    
    fclose(Fid);
    
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
        save(['FixedNoteFiles\', NoteFiles(i).name], 'labels', 'Fs', 'offsets', 'onsets', 'ONMultiplier', 'OFFMultiplier');    
    else
        save(['FixedNoteFiles/', NoteFiles(i).name], 'labels', 'Fs', 'offsets', 'onsets', 'ONMultiplier', 'OFFMultiplier');
    end
end
disp('Finished saving fixed note files');