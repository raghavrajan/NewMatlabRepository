function [] = PlotSyllablesBeforeMotif(NoteFiles, DataDir, FileType, Motif, InterBoutInterval)

for i = 1:length(NoteFiles),
    Notes = load(NoteFiles(i).name);
    if (strfind(FileType,'obs'))
        [Song,Fs] = soundin_copy(DataDir, NoteFiles(i).name(1:end-8), 'obs0r');
        Song = Song * 5/32768;
    else
        if (strfind(FileType,'wav'));
            cd(DataDir);
            [Song, Fs] = wavread(NoteFiles(i).name(1:end-8));
            cd(PresentDir);
        else 
            if (strfind(FileType, 'okrank'))
                [Song, Fs] = ReadOKrankData(DataDir, NoteFiles(i).name(1:end-8), 1);
            end
        end
    end
    
    F_High = 1000; % high pass filter cut-off
    F_Low = 7000; % low pass filter cut-off
    
    FilterForSong = fir1(80, [F_High*2/Fs F_Low*2/Fs], 'bandpass');
    FiltSong = filtfilt(FilterForSong, 1, Song);
    SmoothWinSize = 0.008;
    
    Window = ones(round(SmoothWinSize*Fs), 1);
    Window = Window/sum(Window);
    smooth = 10*log10(conv(FiltSong.*FiltSong, Window, 'same'));
    SongTime = (0:1/Fs:(length(Song)-1)/Fs);
    
    NewFs = 1000;
    Time = (0:1/NewFs:(length(Song)-1)/Fs);
    Time(end) = [];
    smooth = spline(SongTime, smooth, Time);
    smooth = smooth + abs(min(smooth));
    
    InterSyllIntervals = (Notes.onsets(2:end) - Notes.offsets(1:end-1))/1000;
    Bouts = find(InterSyllIntervals > InterBoutInterval);
    
    BoutIndices = [];
    BoutLengths = [];
    for i = 0:length(Bouts),
        if (i ~= length(Bouts))
            if (i == 0)
                BoutIndices = [BoutIndices; [1 Bouts(i+1)]];
                BoutLengths = [BoutLengths; (Notes.offsets(Bouts(i+1)) - Notes.onsets(1))/1000];
            else
                BoutIndices = [BoutIndices; [(Bouts(i) + 1) Bouts(i+1)]];
                BoutLengths = [BoutLengths; (Notes.offsets(Bouts(i+1)) - Notes.onsets(Bouts(i) + 1))/1000];
            end
        else
            if (i == 0)
                BoutIndices = [BoutIndices; [1 length(Notes.onsets)]];
                BoutLengths = [BoutLengths; (Notes.offsets(end) - Notes.onsets(1))/1000];
            else
                BoutIndices = [BoutIndices; [(Bouts(i) + 1) length(Notes.onsets)]];
                BoutLengths = [BoutLengths; (Notes.offsets(end) - Notes.onsets(Bouts(i) + 1))/1000];
            end
        end
    end
    
    for j = 1:size(Bouts,1),
        if (Notes.onsets(Bouts(j,1)
    Matches = [];
    for j = 1:length(Motif),
        Matches = [Matches find(Notes.labels == Motif{j}(1))];
    end
    
    for j = 1:length(Matches),
        if (Notes.onsets(
end