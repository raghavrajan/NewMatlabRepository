function [MotifIntroNoteIndices, NonMotifIntroNoteIndices] = INACalcSAPFeatures(MotifIntroNoteIndices, NonMotifIntroNoteIndices, Onsets, Offsets, SongFile, FileType, DataDir, RecFileDir)

PresentDir = pwd;

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
    
Time = (0:1:(length(Song)-1))/Fs;

[m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude , m_Freq, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight] = deriv(Song, Fs);

% FFTWinSize = 0.002;
% WinSize = round(FFTWinSize * Fs);
% WinOverlap = WinSize - 4;
% [S, F, T, P] = spectrogram(Song, hamming(WinSize), WinOverlap, WinSize, Fs);
% Freq = find((F >= 860) & (F <= 8600));
% smooth = log10(sum(S(Freq,:).*conj(S(Freq,:))));
% Fs1 = 1/(T(2) - T(1));
% Width = 0.004;
% GaussianLen = 1;
% XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs1)));
% XGauss = XGauss - (length(XGauss) + 1)/2;
% GaussWin = (1/((Width * Fs1) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs1) * (Width * Fs1)));
% smooth = conv(smooth, GaussWin, 'same');
% 
% Obj = gmdistribution.fit(smooth', 2, 'Start', 'randSample', 'Replicates', 10);
% if (Obj.mu(1) > Obj.mu(2))
%     Threshold = Obj.mu(1);
% else
%     Threshold = Obj.mu(2);
% end
% 
% Amplitude = spline(T, smooth, Time);
% [onsets, offsets] = segment_song_for_labelling(Amplitude, Fs, 1, 1, Threshold, [Obj.mu(1) Obj.Sigma(1)]);
% onsets = onsets/1000;
% offsets = offsets/1000;

T = linspace(Time(1), Time(end), length(m_Entropy));
if (~isfield(MotifIntroNoteIndices, 'Onsets'))
    return;
end
for i = 1:length(MotifIntroNoteIndices.Onsets),
    for j = 1:length(MotifIntroNoteIndices.Onsets{i}),
%         [Val, Ind] = min(abs(onsets - MotifIntroNoteIndices.Onsets{i}(j)));
%         if (Val <= 0.005)
%             MotifIntroNoteIndices.Onsets{i}(j) = onsets(Ind);
%         end
%         [Val, Ind] = min(abs(offsets - MotifIntroNoteIndices.Offsets{i}(j)));
%         if (Val <= 0.005)
%             MotifIntroNoteIndices.Offsets{i}(j) = offsets(Ind);
%         end
        StartIndex = find(T <= MotifIntroNoteIndices.Onsets{i}(j), 1, 'last');
        EndIndex = find(T >= MotifIntroNoteIndices.Offsets{i}(j), 1, 'first');
        MotifIntroNoteIndices.PG{i}(j) = mean(m_PitchGoodness(StartIndex:EndIndex));
        MotifIntroNoteIndices.Amplitude{i}(j) = mean(m_amplitude(StartIndex:EndIndex));
        MotifIntroNoteIndices.FM{i}(j) = mean(m_FM(StartIndex:EndIndex));
        MotifIntroNoteIndices.AM{i}(j) = mean(m_AM(StartIndex:EndIndex));
        MotifIntroNoteIndices.Pitch{i}(j) = mean(Pitch_chose(StartIndex:EndIndex));
        MotifIntroNoteIndices.Entropy{i}(j) = mean(m_Entropy(StartIndex:EndIndex));
        MotifIntroNoteIndices.MeanFreq{i}(j) = mean(m_Freq(StartIndex:EndIndex));
    end
end

for i = 1:length(MotifIntroNoteIndices.Onsets),
    if (~isempty(MotifIntroNoteIndices.Onsets{i}))
        MotifIntroNoteIndices.Durs{i} = MotifIntroNoteIndices.Offsets{i} - MotifIntroNoteIndices.Onsets{i};
%        [Val, Ind] = min(abs(onsets - Onsets(MotifIntroNoteIndices.Indices{i}(end) + 1)));
        MotifIntroNoteIndices.Intervals{i} = [MotifIntroNoteIndices.Onsets{i}(2:end); Onsets(MotifIntroNoteIndices.Indices{i}(end) + 1)] - MotifIntroNoteIndices.Offsets{i};
    end
end

for i = 1:length(NonMotifIntroNoteIndices.Onsets),
    for j = 1:length(NonMotifIntroNoteIndices.Onsets{i}),
 %       [Val, Ind] = min(abs(onsets - NonMotifIntroNoteIndices.Onsets{i}(j)));
 %       NonMotifIntroNoteIndices.Onsets{i}(j) = onsets(Ind);
 %       [Val, Ind] = min(abs(offsets - NonMotifIntroNoteIndices.Offsets{i}(j)));
 %       NonMotifIntroNoteIndices.Offsets{i}(j) = offsets(Ind);

        StartIndex = find(T <= NonMotifIntroNoteIndices.Onsets{i}(j), 1, 'last');
        EndIndex = find(T >= NonMotifIntroNoteIndices.Offsets{i}(j), 1, 'first');
        NonMotifIntroNoteIndices.PG{i}(j) = mean(m_PitchGoodness(StartIndex:EndIndex));
        NonMotifIntroNoteIndices.Amplitude{i}(j) = mean(m_amplitude(StartIndex:EndIndex));
        NonMotifIntroNoteIndices.FM{i}(j) = mean(m_FM(StartIndex:EndIndex));
        NonMotifIntroNoteIndices.AM{i}(j) = mean(m_AM(StartIndex:EndIndex));
        NonMotifIntroNoteIndices.Pitch{i}(j) = mean(Pitch_chose(StartIndex:EndIndex));
        NonMotifIntroNoteIndices.Entropy{i}(j) = mean(m_Entropy(StartIndex:EndIndex));
        NonMotifIntroNoteIndices.MeanFreq{i}(j) = mean(m_Freq(StartIndex:EndIndex));
    end
end

for i = 1:length(NonMotifIntroNoteIndices.Onsets),
    if (~isempty(NonMotifIntroNoteIndices.Onsets{i}))
        NonMotifIntroNoteIndices.Durs{i} = NonMotifIntroNoteIndices.Offsets{i} - NonMotifIntroNoteIndices.Onsets{i};
    end
end