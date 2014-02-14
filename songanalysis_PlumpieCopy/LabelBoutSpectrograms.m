function [Bout] = LabelBoutSpectrograms(DirectoryName, FileList, FileType, PlotOption, NoiseMSD, SoundMSD, ThreshMultiplier, MotifTemplate, MotifTemplateLabels, StretchValues)

Threshold = (NoiseMSD(1) + SoundMSD(1))/2;
disp(['Threshold is ', num2str(Threshold)]);

min_int = 0.003;
min_dur = 0.001;
sm_win = 0.008;

FFTWinSize = sm_win; % in sec
FFTWinOverlap = 0.9;

if ispc
    if (DirectoryName(end) ~= '\')
        DirectoryName(end + 1) = '\';
    end 
else
    if (DirectoryName(end) ~= '/')
        DirectoryName(end + 1) = '/';
    end
end

DiscardedBouts = 0;

BoutIndex = 1;

Fid = fopen(FileList, 'r');
SongFile = fgetl(Fid);

cd(DirectoryName);

while (ischar(SongFile(1)))
    StartBoutIndex = BoutIndex;
    Slash = find((SongFile == '/') | (SongFile == '\'));
    if (~isempty(Slash))
        SongFile = SongFile(Slash(end)+1:end);
    end
    disp(SongFile);
    try
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
    catch
        continue;
    end
    if (isempty(Song))
        continue;
    end
    
    Time = (1:1:length(Song))/Fs;
    
    FFTWinSize = 0.002;
    FFTWinOverlap = 0.98;
    WinSize = round(FFTWinSize * Fs);
    %WinOverlap = round(FFTWinOverlap * WinSize);
    WinOverlap = WinSize - 4;
    
    [S1, F, T, P] = spectrogram(Song, hamming(WinSize), WinOverlap, WinSize, Fs);
    
    Freq1 = find((F >= 860) & (F <= 7000));
    S1 = log10(abs(S1(Freq1,:)));
    Power = sum(S1);
    Fs1 = 1/(T(2) - T(1));
    
    % Power = 10*log10(sum(P));
    % ThreshCrossing = Power > Threshold;
    Width = 0.005;
    GaussianLen = 4;
    XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs1)));
    XGauss = XGauss - (length(XGauss) + 1)/2;
    GaussWin = (1/((Width * Fs1) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs1) * (Width * Fs1)));
    Power = conv(Power, GaussWin, 'same');
    
    try
        [onsets, offsets] = segment_song_for_labelling(Power, (1/(T(2) - T(1))), min_int*1000, min_dur*1000, Threshold, NoiseMSD);
        onsets = onsets/1000;
        offsets = offsets/1000;
    catch
        disp(['No onsets in ', SongFile]);
        clear onsets offsets Bouts BoutOnsets BoutOffsets Temp;
        SongFile = fgetl(Fid);
        continue;
    end

    if (size(onsets,1) < size(onsets, 2))
        onsets = onsets';
    end
    if (size(offsets,1) < size(offsets, 2))
        offsets = offsets';
    end
    
    Temp = (onsets(2:end) - offsets(1:end-1));
    if (onsets(1) > 0.5)
        BoutOnsets = [onsets(1); onsets(find(Temp > 1) + 1)];
        if (offsets(end) > (T(end) - 0.5))
            BoutOffsets = [offsets(Temp > 1); (offsets(end) - 0.5)];
        else
            BoutOffsets = [offsets(Temp > 1); offsets(end)];
        end
    else
        BoutOnsets = [(onsets(1) + 0.5); onsets(find(Temp > 1) + 1)];
        if (offsets(end) > (T(end) - 0.5))
            BoutOffsets = [offsets(Temp > 1); (offsets(end) - 0.5)];
        else
            BoutOffsets = [offsets(Temp > 1); offsets(end)];
        end
    end
    BoutOnsets = BoutOnsets - 0.5;
    BoutOffsets = BoutOffsets + 0.5;

    AmpDeriv = diff(Power);
    AmpDeriv = (AmpDeriv - mean(AmpDeriv))/std(AmpDeriv);
    AmpDeriv = spline(T(1:end-1), AmpDeriv, Time);
    
    [pks, locs] = findpeaks(AmpDeriv, 'MinPeakHeight', 1.5);
    onsets = Time(locs);
    DerivOnsets = Time(locs);
    
    [pks, locs] = findpeaks(-AmpDeriv, 'MinPeakHeight', 2);
    offsets = Time(locs);
    DerivOffsets = Time(locs);
    
    FFTWinSize = sm_win;
    FFTWinOverlap = 0.5;
    WinSize = round(FFTWinSize * Fs);
    WinOverlap = round(FFTWinOverlap * WinSize);

    [S, F, T1, P] = spectrogram(Song, hamming(WinSize), WinOverlap, WinSize, Fs);
    
    Freq1 = find((F >= 860) & (F <= 7000));
    S = log10(abs(S(Freq1,:)));
    
    Bouts = [BoutOnsets BoutOffsets];
    if (~isempty(Bouts))
        DiscardedBouts = DiscardedBouts + length(find((Bouts(:,2) - Bouts(:,1)) <= 2));
        Bouts = Bouts((Bouts(:,2) - Bouts(:,1)) > 2,:);

        for BoutNo = 1:size(Bouts,1),
            Indices = (find((T1 >= Bouts(BoutNo,1)) & (T1 <= Bouts(BoutNo, 2))));

            TempS = S(:,Indices);
            for Sylls = 1:length(MotifTemplate),
                for OnOffs = 1:length(MotifTemplate{Sylls}),
                    x = 1:1:size(MotifTemplate{Sylls}{OnOffs}, 2);
                    WarpIndex = 1;
                    for i = [StretchValues],
                        xx = linspace(x(1), x(end), round(size(MotifTemplate{Sylls}{OnOffs},2) * (1 + i/100)));
                        clear WMotif;
                        WMotif = zeros(size(MotifTemplate{Sylls}{OnOffs},1), length(xx));
                        for j = 1:size(MotifTemplate{Sylls}{OnOffs}, 1);
                            WMotif(j,:) = spline(x, MotifTemplate{Sylls}{OnOffs}(j,:), xx);
                        end
                        TempS = S(:,Indices);

                        WinMean = zeros((size(TempS,2) - size(WMotif, 2) + 1), 1);
                        WinSTD = zeros((size(TempS,2) - size(WMotif, 2) + 1), 1);

                        for ColNo = 1:(size(TempS,2) - size(WMotif, 2) + 1),
                            StartIndex = ((ColNo - 1)*size(WMotif,1)) + 1;
                            WinIndices = StartIndex:1:(StartIndex + size(WMotif,1)*size(WMotif,2) - 1);
                            WinMean(ColNo) = mean(TempS(WinIndices));
                            WinSTD(ColNo) = std(TempS(WinIndices));
                        end
                        [Match] = CalTemplateMatch(WMotif, TempS, WinMean, WinSTD);
                        %Match = (Match - mean(mean(Match)))/std(Match);                    
                        TempSeqMatch{WarpIndex} = Match;
                        WarpIndex = WarpIndex + 1;
                        fprintf('>');
                    end
                    fprintf('\t');
                    clear Match;
                    for i = 1:length(TempSeqMatch),
                        Match(i,:) = TempSeqMatch{i}(1:length(TempSeqMatch{end}));
                    end
                    Match = max(Match);
                    Match = (Match - mean(mean(Match)))/std(Match);
                    Bout.BoutSeqMatch{BoutIndex}{Sylls, OnOffs} = Match;            
                    clear TempSeqMatch;
                    [pks, SyllOnsets] = findpeaks(Match, 'MinPeakHeight', 4);
                    Bout.BoutSyllOnsets{BoutIndex}{Sylls, OnOffs} = SyllOnsets;
                    Indices2 = (find((T >= Bouts(BoutNo,1)) & (T <= Bouts(BoutNo, 2))));
                    TempT = T(Indices2);
                    TempT2 = T1(Indices);
                    Power1 = sum(S(:,Indices));
                    xx = linspace(TempT(1), TempT(end), length(TempT)*10);
                    Power2 = spline(TempT, Power(Indices2), xx);
                    PlusSD = xx(find(Power2 > (NoiseMSD(1) + ThreshMultiplier*NoiseMSD(2))));
                    MinusSD = xx(find(Power2 <= (NoiseMSD(1) + ThreshMultiplier*NoiseMSD(2))));
                    RealOnset = [];
                    RealLabel = [];
                    RealOffset = [];
                    for SyllOnNo = 1:length(SyllOnsets),
                        if (OnOffs == 1)
                            [minval, minind] = min(abs(onsets - TempT2(SyllOnsets(SyllOnNo))));
                            RealOnset(SyllOnNo) = onsets(minind);
                        else
                            [minval, minind] = min(abs(offsets - TempT2(SyllOnsets(SyllOnNo))));
                            RealOffset(SyllOnNo) = offsets(minind);
                        end
                        RealLabel = [RealLabel MotifTemplateLabels{Sylls}{OnOffs}];
                    end
                    Bout.BoutOnsets{BoutIndex}{Sylls, OnOffs} = RealOnset;
                    Bout.BoutOffsets{BoutIndex}{Sylls, OnOffs} = RealOffset;
                    Bout.BoutLabels{BoutIndex}{Sylls, OnOffs} = RealLabel;
                end
            end
            Bout.T{BoutIndex} = T1(Indices);
            RealOnsets = [];
            RealOnIndices = [];
            RealOffsets = [];
            RealLabels = [];
            for Sylls = 1:size(Bout.BoutOnsets{BoutIndex},1),
                RealOnsets = [RealOnsets unique(Bout.BoutOnsets{BoutIndex}{Sylls, 1})];
                RealOnIndices = [RealOnIndices ones(1,length(unique(Bout.BoutOnsets{BoutIndex}{Sylls, 1})))*Sylls];
            end
            [RealOnsets, RealOnsetsIndices] = sort(RealOnsets);
            RealOnIndices = RealOnIndices(RealOnsetsIndices);
            for RealOns = 1:length(RealOnsets),
                if (RealOnIndices(RealOns) == 1)
                    if (RealOns == length(RealOnsets))
                        RealOffsetIndex = find(Bout.BoutOffsets{BoutIndex}{1, 2} >= RealOnsets(RealOns), 1, 'first');
                        if (length(RealOffsetIndex) == 0)
                            Duration = 1000;
                        else
                            Duration = Bout.BoutOffsets{BoutIndex}{1, 2}(RealOffsetIndex) - RealOnsets(RealOns);
                        end
                        
                        if (Duration < (1.5*(size(MotifTemplate{RealOnIndices(RealOns)}{1},2)*(T1(2) - T1(1)))))
                            RealOffsets = [RealOffsets Bout.BoutOffsets{BoutIndex}{1, 2}(RealOffsetIndex)];
                        else
                            [MinOff, MinOffIndex] = min(abs(DerivOffsets - (RealOnsets(RealOns) + (size(MotifTemplate{RealOnIndices(RealOns)}{1},2)*(T1(2) - T1(1))))));
                            RealOffsets = [RealOffsets DerivOffsets(MinOffIndex)];
                        end
                    else
                        if (RealOnIndices(RealOns + 1) == RealOnIndices(RealOns))
                            RealOffsetIndex = find(Bout.BoutOffsets{BoutIndex}{1, 2} >= RealOnsets(RealOns), 1, 'first');
                            if (length(RealOffsetIndex) == 0)
                                Duration = 1000;
                            else
                                Duration = Bout.BoutOffsets{BoutIndex}{1, 2}(RealOffsetIndex) - RealOnsets(RealOns);
                            end
                        
                            if (Duration < (1.5*(size(MotifTemplate{RealOnIndices(RealOns)}{1},2)*(T1(2) - T1(1)))))
                                RealOffsets = [RealOffsets Bout.BoutOffsets{BoutIndex}{1, 2}(RealOffsetIndex)];
                            else
                                [MinOff, MinOffIndex] = min(abs(DerivOffsets - (RealOnsets(RealOns) + (size(MotifTemplate{RealOnIndices(RealOns)}{1},2)*(T1(2) - T1(1))))));
                                RealOffsets = [RealOffsets DerivOffsets(MinOffIndex)];
                            end
                        else
                            RealOffsetIndex = find(Bout.BoutOffsets{BoutIndex}{1, 3} >= RealOnsets(RealOns), 1, 'first');
                            if (length(RealOffsetIndex) == 0)
                                Duration = 1000;
                            else
                                Duration = Bout.BoutOffsets{BoutIndex}{1, 3}(RealOffsetIndex) - RealOnsets(RealOns);
                            end
                        
                            if (Duration < (1.5*(size(MotifTemplate{RealOnIndices(RealOns)}{1},2)*(T1(2) - T1(1)))))
                                RealOffsets = [RealOffsets Bout.BoutOffsets{BoutIndex}{1, 3}(RealOffsetIndex)];
                            else
                                [MinOff, MinOffIndex] = min(abs(DerivOffsets - (RealOnsets(RealOns) + (size(MotifTemplate{RealOnIndices(RealOns)}{1},2)*(T1(2) - T1(1))))));
                                RealOffsets = [RealOffsets DerivOffsets(MinOffIndex)];
                            end
                        end
                    end
                else
                    if (RealOns == length(RealOnsets))
                        RealOffsetIndex = find(Bout.BoutOffsets{BoutIndex}{RealOnIndices(RealOns), 3} >= RealOnsets(RealOns), 1, 'first');
                        if (length(RealOffsetIndex) == 0)
                            Duration = 1000;
                        else
                            Duration = Bout.BoutOffsets{BoutIndex}{RealOnIndices(RealOns), 3}(RealOffsetIndex) - RealOnsets(RealOns);
                        end
                        
                        if (Duration < (1.5*(size(MotifTemplate{RealOnIndices(RealOns)}{1},2)*(T1(2) - T1(1)))))
                            RealOffsets = [RealOffsets Bout.BoutOffsets{BoutIndex}{RealOnIndices(RealOns), 3}(RealOffsetIndex)];
                        else
                            [MinOff, MinOffIndex] = min(abs(DerivOffsets - (RealOnsets(RealOns) + (size(MotifTemplate{RealOnIndices(RealOns)}{1},2)*(T1(2) - T1(1))))));
                            RealOffsets = [RealOffsets DerivOffsets(MinOffIndex)];
                        end
                    else
                        RealOffsetIndex = find(Bout.BoutOffsets{BoutIndex}{RealOnIndices(RealOns), 2} >= RealOnsets(RealOns), 1, 'first');
                        if (length(RealOffsetIndex) == 0)
                            Duration = 1000;
                        else
                            Duration = Bout.BoutOffsets{BoutIndex}{RealOnIndices(RealOns), 2}(RealOffsetIndex) - RealOnsets(RealOns);
                        end
                        
                        if (Duration < (1.5*(size(MotifTemplate{RealOnIndices(RealOns)}{1},2)*(T1(2) - T1(1)))))
                            RealOffsets = [RealOffsets Bout.BoutOffsets{BoutIndex}{RealOnIndices(RealOns), 2}(RealOffsetIndex)];
                        else
                            [MinOff, MinOffIndex] = min(abs(DerivOffsets - (RealOnsets(RealOns) + (size(MotifTemplate{RealOnIndices(RealOns)}{1},2)*(T1(2) - T1(1))))));
                            RealOffsets = [RealOffsets DerivOffsets(MinOffIndex)];
                        end
                    end
                end
                RealLabels = [RealLabels MotifTemplateLabels{RealOnIndices(RealOns)}{1}];
            end
            if (length(RealOnsets) > 0)
                [UniqueOnsets, UniqueOnIndices] = unique(RealOnsets);
                UniqueOffsets = RealOffsets(UniqueOnIndices);
                UniqueLabels = RealLabels(UniqueOnIndices);
            else
                UniqueOnsets = RealOnsets;
                UniqueOffsets = RealOffsets;
                UniqueLabels = RealLabels;
            end
            if (length(RealOnsets) ~= length(UniqueOnsets))
                DiffOnsets = setdiff((1:1:length(RealOnsets)), UniqueOnIndices);
                for DiffOnNo = 1:length(DiffOnsets),
                    MaxVal = -100;
                    DiffIndices = find(RealOnsets == RealOnsets(DiffOnsets(DiffOnNo)));
                    for TempDiffOns = 1:length(DiffIndices),
                        TempIndex = find(Bout.T{BoutIndex} <= RealOnsets(DiffOnsets(DiffOnNo)), 1, 'last');
                        if (TempIndex < (length(Bout.BoutSeqMatch{BoutIndex}{RealOnIndices(DiffIndices(TempDiffOns)),1}) - 10))
                            if (TempIndex <= 10)
                               if (MaxVal < max(Bout.BoutSeqMatch{BoutIndex}{RealOnIndices(DiffIndices(TempDiffOns)), 1}((1):(TempIndex+10))))
                                    MaxInd = RealOnIndices(DiffIndices(TempDiffOns));
                                    DiffOnIndex = DiffIndices(TempDiffOns); 
                                    MaxVal = max(Bout.BoutSeqMatch{BoutIndex}{RealOnIndices(DiffIndices(TempDiffOns)), 1}((1):(TempIndex+10)));
                                end
                            else
                                if (MaxVal < max(Bout.BoutSeqMatch{BoutIndex}{RealOnIndices(DiffIndices(TempDiffOns)), 1}((TempIndex-10):(TempIndex+10))))
                                    MaxInd = RealOnIndices(DiffIndices(TempDiffOns));
                                    DiffOnIndex = DiffIndices(TempDiffOns); 
                                    MaxVal = max(Bout.BoutSeqMatch{BoutIndex}{RealOnIndices(DiffIndices(TempDiffOns)), 1}((TempIndex-10):(TempIndex+10)));
                                end
                            end
                        else
                            if (MaxVal < max(Bout.BoutSeqMatch{BoutIndex}{RealOnIndices(DiffIndices(TempDiffOns)), 1}((TempIndex-10):end)))
                                MaxInd = RealOnIndices(DiffIndices(TempDiffOns));
                                DiffOnIndex = DiffIndices(TempDiffOns); 
                                MaxVal = max(Bout.BoutSeqMatch{BoutIndex}{RealOnIndices(DiffIndices(TempDiffOns)), 1}((TempIndex-10):end));
                            end
                        end
                    end
                    UniqueOffsets(find(UniqueOnsets == RealOnsets(DiffOnsets(DiffOnNo)))) = RealOffsets(DiffOnIndex);
                    UniqueLabels(find(UniqueOnsets == RealOnsets(DiffOnsets(DiffOnNo)))) = RealLabels(DiffOnIndex);
                end
            end
            
            if (length(UniqueOffsets) ~= length(unique(UniqueOffsets)))
                [UniqueOffsets, UniqueOffIndices] = unique(UniqueOffsets, 'first');
                UniqueOnsets = UniqueOnsets(UniqueOffIndices);
                UniqueLabels = UniqueLabels(UniqueOffIndices);
            end
            
            Bout.UniqueOnsets{BoutIndex} = UniqueOnsets;
            Bout.UniqueOffsets{BoutIndex} = UniqueOffsets;
            Bout.UniqueLabels{BoutIndex} = UniqueLabels;
            Bout.SongFileName{BoutIndex} = SongFile;
            BoutIndex = BoutIndex + 1;
            fprintf('\t');
        end
        
    end
    onsets = [];
    offsets = [];
    labels = [];
    for BoutNo = StartBoutIndex:BoutIndex-1,
        onsets = [onsets Bout.UniqueOnsets{BoutNo}];
        offsets = [offsets Bout.UniqueOffsets{BoutNo}];
        labels = [labels Bout.UniqueLabels{BoutNo}];
    end
    onsets = onsets' * 1000;
    offsets = offsets' * 1000;

    if (strfind(PlotOption, 'on'))
        PlotSpectrogram(DirectoryName, SongFile, FileType);
        hold on;
        plot(T, (Power + (abs(min(Power))))*100);
        for j = 1:length(onsets),
            plot([onsets(j)/1000 onsets(j)/1000 offsets(j)/1000 offsets(j)/1000], [0 8000 8000 0],'r');
            text(onsets(j)/1000, 8200, labels(j));
        end
        
        if (~isempty(Bouts))
            for j = 1:size(Bouts,1),
                plot([Bouts(j,1) Bouts(j,2)], [7500 7500],'b', 'LineWidth', 2);
                %plot(Bout.T{BoutIndex - j}(1:length(Bout.BoutSeqMatch{BoutIndex - j})), (Bout.BoutSeqMatch{BoutIndex-j}/max(Bout.BoutSeqMatch{BoutIndex-j})) *8000, 'r', 'LineWidth', 2);
            end
        end
        zoom xon;
        uiwait(gcf);
    else
        fprintf('\n');
    end
    
    if (ispc)
        save(['C:\Documents and Settings\PlumPie\My Documents\RaghavData\', SongFile, '.not.mat'], 'labels', 'onsets', 'offsets');
    else
        save([SongFile, '.not.mat'], 'labels', 'onsets', 'offsets');
    end
    
    fprintf('%s\n', labels);    
    clear onsets offsets Bouts BoutOnsets BoutOffsets Temp;
    SongFile = fgetl(Fid);
end
disp(['Discarded ', num2str(DiscardedBouts), ' bouts out of a total of ', num2str(BoutIndex - 1), ' bouts']);
