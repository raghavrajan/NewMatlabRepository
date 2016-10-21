function [] = ASSLMakeTemplatesSpectralMatchAnalysis(DirName, FileName, FileType, SongChanNo, XDur, FFTWinSize, FFTWinOverlap, OutputFileDir, OutputFileName, SyllableLabel)

[RawData, Fs] = ASSLGetRawData(DirName, FileName, FileType, SongChanNo);
    
Time = (1:1:length(RawData))/Fs;

Motif = RawData(find((Time >= XDur(1)) & (Time <= XDur(2))));
WinSize = round(FFTWinSize/1000 * Fs);
WinOverlap = round(FFTWinOverlap/1000 * WinSize);
[S1, F, T, P] = spectrogram(Motif, hamming(WinSize), WinOverlap, WinSize, Fs);
Freq1 = find((F >= 860) & (F <= 8600));
Power = log10(sum(S1(Freq1,:).*conj(S1(Freq1,:))));
S = log10(abs(S1(Freq1,:)));
S = (S - mean(reshape(S, 1, size(S,1)*size(S,2))))/std(reshape(S, 1, size(S,1)*size(S,2)));
figure
contourf(S);
colorbar
MotifTemplate = S;
save([OutputFileDir, '/', OutputFileName], 'MotifTemplate', 'DirName', 'FileName', 'FileType', 'SongChanNo', 'XDur', 'SyllableLabel');

disp(['Saved templates to ', OutputFileName, ' in directory ', OutputFileDir]);