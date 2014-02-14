function [] = AutoSegmentMicrolesionSongs(Dir, RawDataDir, Date, Condition, FileType, PlotOption, SAPFeats, SAPLabels)

PresentDir = pwd;

cd(Dir);
MaxMatchFiles = dir('*maxMatch*.txt');
!mkdir AutoSegmentNoteFiles

for i = 1:length(MaxMatchFiles),
    if (~isempty(strfind(MaxMatchFiles(i).name, [Date, '_', Condition])))
        Temp = textread(MaxMatchFiles(i).name, '%s', 'delimiter', '\n');
    end
end

MaxFiles = min(25, length(Temp));

% Calculate SAP Features for all these files
cd(RawDataDir);
for i = 1:MaxFiles,
    SongFile = Temp{i};
    disp(SongFile);
    if (strfind(FileType, 'okrank'))
        [Song, Fs] = ReadOKrankData(RawDataDir, SongFile, 1);
    else
        if (strfind(FileType, 'wav'))
            [Song, Fs] = wavread(SongFile);
        else
            if (strfind(FileType, 'obs'))
                channel_string = strcat('obs',num2str(0),'r');
                [Song, Fs] = soundin_copy([RawDataDir, '/'], SongFile, channel_string);

                % Convert to V - 5V on the data acquisition is 32768
                Song = Song * 5/32768;
            end
        end
    end
    
    Time = (0:1:length(Song)-1)/Fs;
    F_High = 1000; % high pass filter cut-off
    F_Low = 7000; % low pass filter cut-off
    
    FilterForSong = fir1(80, [F_High*2/Fs F_Low*2/Fs]);
    FiltSong = filtfilt(FilterForSong, 1, Song);
    
%     FFTWinSize = 0.0025;
%     WinSize = round(FFTWinSize * Fs);
%     WinOverlap = round(WinSize * 0.9);
% 
%     [S, F, T, P] = spectrogram(Song, hamming(WinSize), WinOverlap, WinSize, Fs);
%     
%     Freq = find((F >= 860) & (F <= 8600));
%     smooth = 10*log10(sum((S(Freq,:).*conj(S(Freq,:)))));
%     Fs1 = 1/(T(2) - T(1));
%     
%     Window = ones(round(0.004*Fs1), 1);
%     Window = Window/sum(Window);
% 
%     smooth = conv(smooth, Window, 'same');
%     smooth = conv(smooth, Window, 'same');
%     
    SmoothWinSize = 0.005;
    
    Window = ones(round(SmoothWinSize*Fs), 1);
    Window = Window/sum(Window);
    smooth = 10*log10(conv(FiltSong.^2, Window, 'same'));
    
    Obj = gmdistribution.fit(smooth, 2);
    
    [Val, Index] = min(Obj.mu);
    %UpperThreshold = Obj.mu(Index) + 4*sqrt(Obj.Sigma(Index));
    UpperThreshold = mean(Obj.mu);
    %RestValues = smooth(find(smooth < UpperThreshold));
    %LowerThreshold = mean(RestValues) + 2*std(RestValues);
%    LowerThreshold = Obj.mu(Index) - 2*sqrt(Obj.Sigma(Index));

    TempCrossings = zeros(size(Time));
    TempCrossings(find(smooth >= UpperThreshold)) = 1;
    Win = [1 -1];
    Trans = conv(TempCrossings, Win, 'same');
    UpperThreshOnsets = find(Trans>0);
    UpperThreshOffsets = find(Trans<0);
    
%     TempCrossings = zeros(size(Time));
%     TempCrossings(find(smooth >= LowerThreshold)) = 1;
%     Win = [1 -1];
%     Trans = conv(TempCrossings, Win, 'same');
%     LowerThreshOnsets = find(Trans>0);
%     LowerThreshOffsets = find(Trans<0);
   
    Onsets = [];
    Offsets = [];
    OnsetIndex = 1;
    StartTime = 0;
    if (length(UpperThreshOffsets) > length(UpperThreshOnsets))
        UpperThreshOffsets(1) = [];
    else
        if (length(UpperThreshOffsets) < length(UpperThreshOnsets))
            UpperThresOnsets(end) = [];
        end
    end
    
%     for j = 1:length(LowerThreshOnsets),
%         if (LowerThreshOnsets(j) < StartTime)
%             continue;
%         end
%         Diff1 = UpperThreshOnsets(find(UpperThreshOnsets > LowerThreshOnsets(j))) - LowerThreshOnsets(j);
%         Diff2 = LowerThreshOnsets(find(LowerThreshOnsets > LowerThreshOnsets(j))) - LowerThreshOnsets(j);
%         if (~isempty(Diff2) && ~isempty(Diff1))
%             if (min(Diff1) < min(Diff2))
% %               Onsets(OnsetIndex) = Time(LowerThreshOnsets(j));
%                 Onsets(OnsetIndex) = Time(UpperThreshOnsets(find(UpperThreshOnsets > LowerThreshOnsets(j), 1, 'first')));
%                 TempIndex = LowerThreshOffsets(find(LowerThreshOffsets > LowerThreshOnsets(j), 1, 'first'));
%                 Offsets(OnsetIndex) = Time(UpperThreshOffsets(find(UpperThreshOffsets < TempIndex, 1, 'last')));
%                 %Offsets(OnsetIndex) = Time(LowerThreshOffsets(find(LowerThreshOffsets > LowerThreshOnsets(j), 1, 'first')));
%                 OnsetIndex = OnsetIndex + 1;
%                 StartTime = Offsets(OnsetIndex-1);
%             end
%         else
%             if (isempty(Diff2))
%                 Onsets(OnsetIndex) = Time(LowerThreshOnsets(j));
%                 Onsets(OnsetIndex) = Time(UpperThreshOnsets(find(UpperThreshOnsets > LowerThreshOnsets(j), 1, 'first')));
%                 TempIndex = LowerThreshOffsets(find(LowerThreshOffsets > LowerThreshOnsets(j), 1, 'first'));
%                 Offsets(OnsetIndex) = Time(UpperThreshOffsets(find(UpperThreshOffsets < TempIndex, 1, 'last')));
%                 %Offsets(OnsetIndex) = Time(LowerThreshOffsets(find(LowerThreshOffsets > LowerThreshOnsets(j), 1, 'first')));
%                 OnsetIndex = OnsetIndex + 1;
%                 StartTime = Offsets(OnsetIndex-1);
%             end
%         end
%     end
    Onsets = Time(UpperThreshOnsets);
    Offsets = Time(UpperThreshOffsets);
    SmallDurs = find((Offsets - Onsets) <= 0.007);
    Onsets(SmallDurs) = [];
    Offsets(SmallDurs) = [];
    SmallInts = find((Onsets(2:end) - Offsets(1:end-1)) <= 0.003);
    Onsets(SmallInts+1) = [];
    Offsets(SmallInts) = [];
    onsets = Onsets*1000;
    offsets = Offsets*1000;
    [Feats] = CalculateSAPFeatsWithOnsets(Song, Time, Fs, Onsets, Offsets);
    
    if (strfind(PlotOption, 'on'))
        PlotSpectrogram([RawDataDir, '/'], SongFile, FileType)
        for j = 1:length(Onsets),
            plot([Onsets(j) Onsets(j) Offsets(j) Offsets(j)], [0 7000 7000 0], 'b', 'LineWidth', 2);
        end
        uiwait(gcf);
    end    
    labels = repmat('0', 1, length(onsets));
    save([Dir, '/AutoSegmentNoteFiles/', SongFile, '.not.mat'], 'onsets', 'offsets', 'labels', 'UpperThreshold', 'F_High', 'F_Low', 'SmoothWinSize');
end
cd(PresentDir);
