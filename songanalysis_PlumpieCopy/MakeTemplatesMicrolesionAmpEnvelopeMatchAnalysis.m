function [] = MakeTemplatesMicrolesionAmpEnvelopeMatchAnalysis(RawSong, Fs, XDur, OutputFileDir, OutputFileName)

Time = (0:1:length(RawSong)-1)/Fs;
sm_win = 0.008;
FFTWinSize = sm_win; % in sec
FFTWinOverlap = 0.5;
Motif = RawSong(find((Time >= XDur(1)) & (Time <= XDur(2))));
WinSize = round(FFTWinSize * Fs);
WinOverlap = round(FFTWinOverlap * WinSize);
[S1, F, T, P] = spectrogram(Motif, hamming(WinSize), WinOverlap, WinSize, Fs);
Freq1 = find((F >= 860) & (F <= 8600));
Power = log10(sum(S1(Freq1,:).*conj(S1(Freq1,:))));
S = sum(log10(abs(S1(Freq1,:))));
S = (S - mean(S))/std(S);
figure
plot(S);
MotifTemplate = S;
save([OutputFileDir, '/', OutputFileName], 'MotifTemplate');
disp(['Saved templates to ', OutputFileName, ' in directory ', OutputFileDir]);