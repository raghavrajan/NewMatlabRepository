function [Motif, Labels] = MakeTemplatesForLabelling(DirectoryName, SongFile, FileType)

cd(DirectoryName);

disp(SongFile);
if (strfind(FileType, 'okrank'))
    [Song, Fs] = ReadOKrankData(DirectoryName, SongFile, 1);
else
    if (strfind(FileType, 'wav'))
        [Song, Fs] = wavread(SongFile);
    else
        if (strfind(FileType, 'obs'))
            channel_string = strcat('obs',num2str(0),'r');
            [Song, Fs] = soundin_copy(DirectoryName, SongFile, channel_string);

            % Convert to uV - 5V on the data acquisition is 32768
            Song = Song * 5/32768;
        end
    end
end

FFTWinSize = 0.008;
FFTWinOverlap = 0.5;
WinSize = round(FFTWinSize * Fs);
WinOverlap = round(FFTWinOverlap * WinSize);

Notes = load([SongFile, '.not.mat']);

Notes.onsets = Notes.onsets/1000;
Notes.offsets = Notes.offsets/1000;

Time = (1:1:length(Song))/Fs;

UniqueSylls = unique(Notes.labels);
for i = 1:length(UniqueSylls),
    Matches = find(Notes.labels == UniqueSylls(i));
    if (UniqueSylls(i) == 'i');
        Matches = Matches(end);
    else
        Matches = Matches(1);
    end
    Syll = Song(find((Time >= (Notes.onsets(Matches) - 0.01)) & (Time <= (Notes.offsets(Matches) + 0.01))));
    [S, F, T1, P] = spectrogram(Syll, hamming(WinSize), WinOverlap, WinSize, Fs);
    Freq1 = find((F >= 860) & (F <= 7000));
    S = log10(abs(S(Freq1,:)));
    S = (S - mean(reshape(S, 1, (size(S,1)*size(S,2)))))/std(reshape(S, 1, (size(S,1)*size(S,2))));
    figure;
    contourf(S);
    title(UniqueSylls(i));
    colorbar;
    Motif{i}{1} = S;
    Labels{i}{1} = UniqueSylls(i);
    
    Syll = Song(find((Time >= (Notes.offsets(Matches) - (Notes.offsets(Matches) - Notes.onsets(Matches))/4)) & (Time <= (Notes.onsets(Matches + 1) + (Notes.offsets(Matches+1) - Notes.onsets(Matches+1))/4)))); 
    [S, F, T1, P] = spectrogram(Syll, hamming(WinSize), WinOverlap, WinSize, Fs);
    Freq1 = find((F >= 860) & (F <= 7000));
    S = log10(abs(S(Freq1,:)));
    S = (S - mean(reshape(S, 1, (size(S,1)*size(S,2)))))/std(reshape(S, 1, (size(S,1)*size(S,2))));
    figure;
    contourf(S);
    colorbar;
    title([UniqueSylls(i), Notes.labels(Matches+1)]);
    Motif{i}{2} = S;
    Labels{i}{2} = [UniqueSylls(i), Notes.labels(Matches+1)];
    
    Syll = Song(find((Time >= (Notes.offsets(Matches) - (Notes.offsets(Matches) - Notes.onsets(Matches))/4)) & (Time <= (Notes.onsets(Matches + 1))))); 
    [S, F, T1, P] = spectrogram(Syll, hamming(WinSize), WinOverlap, WinSize, Fs);
    Freq1 = find((F >= 860) & (F <= 7000));
    S = log10(abs(S(Freq1,:)));
    S = (S - mean(reshape(S, 1, (size(S,1)*size(S,2)))))/std(reshape(S, 1, (size(S,1)*size(S,2))));
    figure;
    contourf(S);
    colorbar;
    title([UniqueSylls(i), ' ']);
    Motif{i}{3} = S;
    Labels{i}{3} = [UniqueSylls(i), ' '];
end