function [] = IntroNoteAmplitudePlots(INoteFile, AlignmentPoint, PreTime, PostTime, ColorScale)

% Alignment Point is a number specifying the point where different bouts
% have to be aligned.
% 0 - refers to alignment on the onset of the first syllable of the motif
% -1 - refers to alignment on the onset of the last intro note
% -2 - refers to alignment on the onset of the second last intro
% -1000 - refers to alignment on the first intro note - in this case no
% trials will be dropped
% note
% and so on .....
%
% If the trial does not have the required number of intro notes for
% alignment, it will be excluded. For eg. if the alignment point is set at
% 'secondlast' and a bout has only one intro note then it will be excluded

% ======================================================================= %

load(INoteFile);

PresentDir = pwd;

TotalBouts = 0;
IFRdt = 0.001;

Width = 0.005;
GaussianLen = 4;
Fs = 1000;

XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * Fs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs) * (Width * Fs)));

if (ispc)
    SlashIndex = find(DataDir == '/');
    ExtIndex = strfind(DataDir, 'raghav');
    DataDir(SlashIndex) = '\';
    DataDir = ['F:\RaghavData\', DataDir((ExtIndex+7):end)];
    if (DataDir(end) ~= '\')
        DataDir(end+1) = '\';
    end
else
    if (DataDir(end) ~= '/')
        DataDir(end+1) = '/';
    end
end

for i = 1:length(MotifIntroNoteIndices),
%    Notes = load([NoteFileDir, '/', ANotes.FileName{i}, '.not.mat']);
    
    if (strfind(FileType,'obs'))
        [Song,Fs] = soundin_copy(DataDir, ANotes.FileName{i}, 'obs0r');
        Song = Song * 5/32768;
    else
        if (strfind(FileType,'wav'));
            cd(DataDir);
            [Song, Fs] = wavread(ANotes.FileName{i});
            cd(PresentDir);
        else 
            if (strfind(FileType, 'okrank'))
                [Song, Fs] = ReadOKrankData(DataDir, ANotes.FileName{i}, 1);
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
    SongTime = (0:1/Fs:ANotes.FileLength{i});
    
    NewFs = 1000;
    Time = (0:1/NewFs:ANotes.FileLength{i});
    Time(end) = [];
    smooth = spline(SongTime, smooth, Time);
    smooth = smooth + abs(min(smooth));
%     FFTWinSize = 0.008;
%     WinSize = round(FFTWinSize * Fs);
%     WinOverlap = round(WinSize * 0.5);
% 
%     [S, F, T, P] = spectrogram(Song, hamming(WinSize), WinOverlap, WinSize, Fs);
%     
%     Freq = find((F >= 860) & (F <= 8600));
%     smooth = 10*log10(sum((S(Freq,:).*conj(S(Freq,:)))));
%     
    INLabels = find(ANotes.Labels{i} == 'i');
    
    for j = INLabels,
        smooth(ceil(ANotes.Onsets{i}(j)*NewFs):floor(ANotes.Offsets{i}(j)*NewFs)) = smooth(ceil(ANotes.Onsets{i}(j)*NewFs):floor(ANotes.Offsets{i}(j)*NewFs))*1.5;
    end

    if (~isempty(MotifIntroNoteIndices{i}))
        UniqueBouts = unique(cell2mat(MotifIntroNoteIndices{i}.BoutNo));
        for k = 1:length(UniqueBouts),
            BoutIndex = find(cell2mat(MotifIntroNoteIndices{i}.BoutNo) == UniqueBouts(k));
            BoutIndex = BoutIndex(1);
            if (MotifIntroNoteIndices{i}.NoofINs{BoutIndex} == 0)
                continue;
            end
            if (MotifIntroNoteIndices{i}.Indices{BoutIndex}(1) == BoutIndices{i}((UniqueBouts(k)-1)*2 + 1))
                if ((AlignmentPoint > -1000) && (MotifIntroNoteIndices{i}.NoofINs{BoutIndex} < abs(AlignmentPoint)))
                    continue;
                end
                TotalBouts = TotalBouts + 1;
                MotifOnset = MotifIntroNoteIndices{i}.MotifOnsets{BoutIndex};
                if (AlignmentPoint == 0)
                    AlignmentTime = MotifOnset;
                else
                    if (AlignmentPoint > -1000)
                        AlignmentTime = MotifIntroNoteIndices{i}.Onsets{BoutIndex}(length(MotifIntroNoteIndices{i}.Onsets{BoutIndex}) - (abs(AlignmentPoint)) + 1);
                    else
                        AlignmentTime = MotifIntroNoteIndices{i}.Onsets{BoutIndex}(1);
                    end
                end
                if (AlignmentTime < (ANotes.FileLength{i} + PostTime))
                    EndTime = AlignmentTime + PostTime;
                else
                    EndTime = ANotes.FileLength{i};
                end

                if (AlignmentTime > PreTime)
                    StartTime = AlignmentTime - PreTime;
                else
                    StartTime = 0;
                end

                TempAmplitude = ones(size(Time));
                
                for j = 1:length(ANotes.Onsets{i}),
                    if (ANotes.Onsets{i}(j) == 0)
                        TempAmplitude(ceil(0.001 * NewFs):floor(ANotes.Offsets{i}(j) * NewFs)) = 0;
                    else
                        TempAmplitude(ceil(ANotes.Onsets{i}(j) * NewFs):floor(ANotes.Offsets{i}(j) * NewFs)) = 0;
                    end
                end
                
                Amplitudes{TotalBouts} = [(Time((find(Time > StartTime, 1, 'first')):(find(Time < EndTime, 1, 'last'))) - AlignmentTime); smooth((find(Time > StartTime, 1, 'first')):(find(Time < EndTime, 1, 'last')))];
                % Amplitudes{TotalBouts} = [(Time((find(Time > StartTime, 1, 'first')):(find(Time < EndTime, 1, 'last'))) - AlignmentTime); TempAmplitude((find(Time > StartTime, 1, 'first')):(find(Time < EndTime, 1, 'last')))];
                NoofINs(TotalBouts) = MotifIntroNoteIndices{i}.NoofINs{BoutIndex};
            end
        end
    end
end

Lengths = cellfun(@length, Amplitudes);
AllAmps = ones(length(Amplitudes), max(Lengths));

% AmpTime = 1:-1/(1/(T(2)-T(1))):(-max(Lengths)/(1/(T(2) - T(1))) + 1);
AmpTime = -PreTime:1/NewFs:PostTime;
AmpTime(end) = [];

[SortedINs, SortedIndices] = sort(NoofINs, 'descend');

for i = 1:length(SortedIndices),
    AllAmps(i,:) = interp1(Amplitudes{SortedIndices(i)}(1,:), Amplitudes{SortedIndices(i)}(2,:), AmpTime);
end
figure
imagesc(AmpTime, [1:1:size(AllAmps,1)], imcomplement(AllAmps/max(max(AllAmps))));
set(gcf, 'Color', 'k');
colormap(ColorScale);
set(gca, 'XColor', 'w', 'YColor', 'w');
xlabel('Time from the beginning of motif (sec)', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'w');
ylabel('Bout #', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'w');
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'w');
colorbar;
disp('Finished Analysis');
