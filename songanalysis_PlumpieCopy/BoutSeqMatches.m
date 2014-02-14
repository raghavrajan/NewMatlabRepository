function [Bout] = BoutSeqMatches(DirectoryName, FileList, FileType, PlotOption, ChooseBouts, min_int, min_dur, GetRidNoise, Template, StretchValues, MinPeakHeight)

sm_win = 0.008;

FFTWinSize = sm_win; % in sec
FFTWinOverlap = 0.9;

if (DirectoryName(end) ~= '/')
    DirectoryName(end + 1) = '/';
end

DiscardedBouts = 0;

BoutIndex = 1;

Fid = fopen(FileList, 'r');
SongFile = fgetl(Fid);

cd(DirectoryName);

while (ischar(SongFile(1)))
    Slash = find((SongFile == '\') | (SongFile == '/'));
    
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
    try
        SongTime = (1:1:length(Song))/Fs;
        [m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude , m_Freq, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight ]=deriv(Song, Fs);
        
        EntropyThreshold = gmdistribution.fit(m_Entropy, 2);
        
        FFTWinSize = sm_win; % in sec
        FFTWinOverlap = 0.9;

        WinSize = round(FFTWinSize * Fs);
        WinOverlap = round(FFTWinOverlap * WinSize);

        [S, F, T, P] = spectrogram(Song, hamming(WinSize), WinOverlap, WinSize, Fs);

        Freq1 = find((F >= 860) & (F <= 8600));
        Power = log10(sum(S(Freq1,:).*conj(S(Freq1,:))));
        % Power = 10*log10(sum(P));
        Power = spline(T, Power, SongTime);

        Obj = gmdistribution.fit(Power', 2);
        [SoundMSD(1), Index] = max(Obj.mu);
        SoundMSD(2) = Obj.Sigma(Index);
        [NoiseMSD(1), Index] = min(Obj.mu);
        NoiseMSD(2) = Obj.Sigma(Index);

        Threshold = (NoiseMSD(1) + SoundMSD(1))/2;
        disp(['Threshold is ', num2str(Threshold)]);

        ThreshCrossing = Power > Threshold;

        [onsets, offsets] = segment_song_for_labelling(Power, Fs, min_int*1000, min_dur*1000, Threshold, NoiseMSD);
        onsets = onsets/1000;
        offsets = offsets/1000;

        if (size(onsets,1) < size(onsets, 2))
            onsets = onsets';
        end
        if (size(offsets,1) < size(offsets, 2))
            offsets = offsets';
        end

        % Get rid of noise based on power in freq. lower than 300Hz
        if (strfind(GetRidNoise, 'yes'))

            NoiseFreq = find(F < 300);
            NoisePower = sum(log10(abs(S(NoiseFreq,:))));
            Width = 0.02;
            XGauss = 1:1:(1 + round(2 * 4 * Width * (1/(T(2) - T(1)))));
            XGauss = XGauss - (length(XGauss) + 1)/2;
            GaussWin = (1/((Width * (1/(T(2) - T(1)))) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * (1/(T(2) - T(1)))) * (Width * (1/(T(2) - T(1))))));
            NTC = NoisePower > 1.5;
            Temp = zeros(size(NTC));
            Temp(find(NTC > 0)) = 1;
            Noise = conv(Temp, GaussWin, 'same');
            Temp2 = zeros(size(Noise));
            Temp2(find(Noise > 0)) = 1;
            Trans = conv(Temp2, [1 -1]);
            NoiseOnsets = T(find(Trans > 0));
            NoiseOffsets = T(find(Trans < 0));

            RemIndices = [];
            for NoiseNo = 1:length(NoiseOnsets),
                RemIndices = [RemIndices; find((onsets >= NoiseOnsets(NoiseNo)) & (onsets <= NoiseOffsets(NoiseNo)))];
                RemIndices = [RemIndices; find((offsets >= NoiseOnsets(NoiseNo)) & (offsets <= NoiseOffsets(NoiseNo)))];
                for (OnsetNo = 1:length(onsets)),
                    if ((onsets(OnsetNo) <= NoiseOnsets(NoiseNo)) && (offsets(OnsetNo) >= NoiseOnsets(NoiseNo)))
                        RemIndices = [RemIndices; OnsetNo];
                    end
                end
            end

            RemIndices = unique(RemIndices);
            onsets(RemIndices) = [];
            offsets(RemIndices) = [];

            disp(['Removed ', num2str(length(RemIndices)), ' syllables that were noise']);
        end
        
        Time = linspace(T(1), T(end), length(m_AM));
        NonEntropicSylls = [];
        for NoteIndex = 1:length(onsets),
            [minval, StartI] = min(abs(Time - onsets(NoteIndex)));
            [minval, EndI] = min(abs(Time - offsets(NoteIndex)));
            if ((EndI - StartI) <= 0)
                continue;
            end
            if (length(find(m_Entropy(StartI:EndI) <= (mean(EntropyThreshold.mu)*2/3))) < (EndI - StartI)/3)
                NonEntropicSylls = [NonEntropicSylls NoteIndex];
            end
        end
        if (~isempty(NonEntropicSylls))
            onsets(NonEntropicSylls) = [];
            offsets(NonEntropicSylls) = [];
        end

        if (strfind(ChooseBouts, 'no'))

            Temp = (onsets(2:end) - offsets(1:end-1));
            if (onsets(1) > 0.5)
                BoutOnsets = [onsets(1); onsets(find(Temp > 1) + 1)];
                if (offsets(end) > (T(end) - 0.5))
                    BoutOffsets = [offsets(find(Temp > 1)); (offsets(end) - 0.5)];
                else
                    BoutOffsets = [offsets(find(Temp > 1)); offsets(end)];
                end
            else
                BoutOnsets = [(onsets(1) + 0.5); onsets(find(Temp > 1) + 1)];
                if (offsets(end) > (T(end) - 0.5))
                    BoutOffsets = [offsets(find(Temp > 1)); (offsets(end) - 0.5)];
                else
                    BoutOffsets = [offsets(find(Temp > 1)); offsets(end)];
                end
            end
            BoutOnsets = BoutOnsets - 0.5;
            BoutOffsets = BoutOffsets + 0.5;

            Bouts = [BoutOnsets BoutOffsets];
            
            FFTWinSize = sm_win;
            FFTWinOverlap = 0.5;
            WinSize = round(FFTWinSize * Fs);
            WinOverlap = round(FFTWinOverlap * WinSize);

            [S, F, T1, P] = spectrogram(Song, hamming(WinSize), WinOverlap, WinSize, Fs);
    
            Freq1 = find((F >= 860) & (F <= 7000));
            S = log10(abs(S(Freq1,:)));
            
            if (length(Bouts) > 0)
                DiscardedBouts = DiscardedBouts + length(find((Bouts(:,2) - Bouts(:,1)) <= 2));
                Bouts = Bouts(find((Bouts(:,2) - Bouts(:,1)) > 2),:);

                Time = linspace(T(1), T(end), length(m_AM));
                for BoutNo = 1:size(Bouts,1),
                    fprintf('B%i> ', BoutNo);
                    OnsetIndices = find((onsets >= Bouts(BoutNo, 1)) & (onsets <= Bouts(BoutNo, 2)));
                    Time = linspace(T(1), T(end), length(m_AM));

                    Indices = (find((T1 >= Bouts(BoutNo,1)) & (T1 <= Bouts(BoutNo, 2))));

                    x = 1:1:size(Template, 2);
                    WarpIndex = 1;
                    for Stretch = [StretchValues],
                        xx = linspace(x(1), x(end), round(size(Template,2) * (1 + Stretch/100)));
                        clear WMotif;
                        WMotif = zeros(size(Template,1), length(xx));
                        for TemplateRow = 1:size(Template, 1);
                            WMotif(TemplateRow,:) = spline(x, Template(TemplateRow,:), xx);
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
                        TempSeqMatch{WarpIndex} = Match;
                        WarpIndex = WarpIndex + 1;
                        fprintf('>');
                    end
                    fprintf('\t');
                    clear Match;
                    for MatchNo = 1:length(TempSeqMatch),
                        Match(MatchNo,:) = TempSeqMatch{MatchNo}(1:length(TempSeqMatch{end}));
                    end
                    Match = max(Match);
                    %Match = (Match - mean(mean(Match)))/std(Match);
                    [pks, SyllMatches] = findpeaks(Match, 'MinPeakHeight', MinPeakHeight);
                    SkippedSylls = 0;
                    if (length(SyllMatches) > 0)
                        for SyllMatchNo = 1:length(SyllMatches),
                            [MinVal, SyllOnsetTime] = min(abs(T1(Indices(SyllMatches(SyllMatchNo))) - onsets(OnsetIndices)));
                            if ((SyllOnsetTime < 3) || (SyllOnsetTime > (length(OnsetIndices) - 2)))
                                SkippedSylls = SkippedSylls + 1;
                            else
                                Bout.Syll.FileName{BoutIndex}{SyllMatchNo} = SongFile;
                                for SyllOrderNo = (SyllOnsetTime - 2):(SyllOnsetTime + 2),
                                    [minval, StartI] = min(abs(Time - onsets(OnsetIndices(SyllOrderNo))));
                                    [minval, EndI] = min(abs(Time - offsets(OnsetIndices(SyllOrderNo))));
                                    Bout.Syll.AM{BoutIndex}(SyllMatchNo, SyllOrderNo + 3 - SyllOnsetTime) = mean(m_AM(StartI:EndI));
                                    Bout.Syll.FM{BoutIndex}(SyllMatchNo, SyllOrderNo + 3 - SyllOnsetTime) = mean(m_FM(StartI:EndI));
                                    Bout.Syll.Entropy{BoutIndex}(SyllMatchNo, SyllOrderNo + 3 - SyllOnsetTime) = mean(m_Entropy(StartI:EndI));
                                    Bout.Syll.Amp{BoutIndex}(SyllMatchNo, SyllOrderNo + 3 - SyllOnsetTime) = mean(m_amplitude(StartI:EndI));
                                    Bout.Syll.Freq{BoutIndex}(SyllMatchNo, SyllOrderNo + 3 - SyllOnsetTime) = mean(m_Freq(StartI:EndI));
                                    Bout.Syll.PG{BoutIndex}(SyllMatchNo, SyllOrderNo + 3 - SyllOnsetTime) = mean(m_PitchGoodness(StartI:EndI));
                                    Bout.Syll.Pitch{BoutIndex}(SyllMatchNo, SyllOrderNo + 3 - SyllOnsetTime) = mean(Pitch_chose(StartI:EndI));
                                    Bout.Syll.Duration{BoutIndex}(SyllMatchNo, SyllOrderNo + 3 - SyllOnsetTime) = Time(EndI) - Time(StartI);
                                end
                            end
                        end
%                             if (SyllOnsetTime == 1)
%                                 Bout.PrevSyll.FileName{BoutIndex}{SyllMatchNo} = SongFile;
%                                 Bout.PrevSyll.AM{BoutIndex}{SyllMatchNo} = 0;
%                                 Bout.PrevSyll.FM{BoutIndex}{SyllMatchNo} = 0;
%                                 Bout.PrevSyll.Entropy{BoutIndex}{SyllMatchNo} = 0;
%                                 Bout.PrevSyll.Amp{BoutIndex}{SyllMatchNo} = 0;
%                                 Bout.PrevSyll.Freq{BoutIndex}{SyllMatchNo} = 0;
%                                 Bout.PrevSyll.PG{BoutIndex}{SyllMatchNo} = 0;
%                                 Bout.PrevSyll.Pitch{BoutIndex}{SyllMatchNo} = 0;
%                                 Bout.PrevSyll.Duration{BoutIndex}{SyllMatchNo} = 0;
%                                 Bout.PrevSyll.GapDur{BoutIndex}{SyllMatchNo} = 0;
%                             else
%                                 [minval, StartI] = min(abs(Time - onsets(OnsetIndices(SyllOnsetTime - 1))));
%                                 [minval, EndI] = min(abs(Time - offsets(OnsetIndices(SyllOnsetTime - 1))));
%                                 Bout.PrevSyll.FileName{BoutIndex}{SyllMatchNo} = SongFile;
%                                 Bout.PrevSyll.AM{BoutIndex}{SyllMatchNo} = mean(m_AM(StartI:EndI));
%                                 Bout.PrevSyll.FM{BoutIndex}{SyllMatchNo} = mean(m_FM(StartI:EndI));
%                                 Bout.PrevSyll.Entropy{BoutIndex}{SyllMatchNo} = mean(m_Entropy(StartI:EndI));
%                                 Bout.PrevSyll.Amp{BoutIndex}{SyllMatchNo} = mean(m_amplitude(StartI:EndI));
%                                 Bout.PrevSyll.Freq{BoutIndex}{SyllMatchNo} = mean(m_Freq(StartI:EndI));
%                                 Bout.PrevSyll.PG{BoutIndex}{SyllMatchNo} = mean(m_PitchGoodness(StartI:EndI));
%                                 Bout.PrevSyll.Pitch{BoutIndex}{SyllMatchNo} = mean(Pitch_chose(StartI:EndI));
%                                 Bout.PrevSyll.Duration{BoutIndex}{SyllMatchNo} = Time(EndI) - Time(StartI);
%                                 Bout.PrevSyll.GapDur{BoutIndex}{SyllMatchNo} = onsets(OnsetIndices(SyllOnsetTime)) - offsets(OnsetIndices(SyllOnsetTime - 1));
%                             end
%                             if (SyllOnsetTime == length(OnsetIndices))
%                                 Bout.NextSyll.FileName{BoutIndex}{SyllMatchNo} = SongFile;
%                                 Bout.NextSyll.AM{BoutIndex}{SyllMatchNo} = 0;
%                                 Bout.NextSyll.FM{BoutIndex}{SyllMatchNo} = 0;
%                                 Bout.NextSyll.Entropy{BoutIndex}{SyllMatchNo} = 0;
%                                 Bout.NextSyll.Amp{BoutIndex}{SyllMatchNo} = 0;
%                                 Bout.NextSyll.Freq{BoutIndex}{SyllMatchNo} = 0;
%                                 Bout.NextSyll.PG{BoutIndex}{SyllMatchNo} = 0;
%                                 Bout.NextSyll.Pitch{BoutIndex}{SyllMatchNo} = 0;
%                                 Bout.NextSyll.Duration{BoutIndex}{SyllMatchNo} = 0;
%                                 Bout.NextSyll.GapDur{BoutIndex}{SyllMatchNo} = 0;
%                             else
%                                 [minval, StartI] = min(abs(Time - onsets(OnsetIndices(SyllOnsetTime + 1))));
%                                 [minval, EndI] = min(abs(Time - offsets(OnsetIndices(SyllOnsetTime + 1))));
%                                 Bout.NextSyll.FileName{BoutIndex}{SyllMatchNo} = SongFile;
%                                 Bout.NextSyll.AM{BoutIndex}{SyllMatchNo} = mean(m_AM(StartI:EndI));
%                                 Bout.NextSyll.FM{BoutIndex}{SyllMatchNo} = mean(m_FM(StartI:EndI));
%                                 Bout.NextSyll.Entropy{BoutIndex}{SyllMatchNo} = mean(m_Entropy(StartI:EndI));
%                                 Bout.NextSyll.Amp{BoutIndex}{SyllMatchNo} = mean(m_amplitude(StartI:EndI));
%                                 Bout.NextSyll.Freq{BoutIndex}{SyllMatchNo} = mean(m_Freq(StartI:EndI));
%                                 Bout.NextSyll.PG{BoutIndex}{SyllMatchNo} = mean(m_PitchGoodness(StartI:EndI));
%                                 Bout.NextSyll.Pitch{BoutIndex}{SyllMatchNo} = mean(Pitch_chose(StartI:EndI));
%                                 Bout.NextSyll.Duration{BoutIndex}{SyllMatchNo} = Time(EndI) - Time(StartI);
%                                 Bout.NextSyll.GapDur{BoutIndex}{SyllMatchNo} = onsets(OnsetIndices(SyllOnsetTime + 1)) - offsets(OnsetIndices(SyllOnsetTime));
%                             end
                    end
                    BoutIndex = BoutIndex + 1;
                    fprintf('\t');
                end
            end
        else
            PlotSpectrogram(DirectoryName, SongFile, FileType);
            zoom xon;
            Flag = 1;
            Bouts = [];
            while(Flag)
                TempMsgBox = msgbox('Zoom into the desired region and then close this box when ready');
                waitfor(TempMsgBox);
                disp('Choose bout onset and offset');
                [x1, y1, button] = ginput(2);
                Bouts(end+1,:) = [x1(1) x1(2)];
                YesNo = inputdlg('Do you want to choose more bouts? (y for yes and n for no)');
                if (isempty(strfind(YesNo{1}, 'y')))
                    Flag = 0;
                end
            end

            Time = linspace(T(1), T(end), length(m_AM));
            for BoutNo = 1:size(Bouts,1),
                fprintf('B%i ', BoutNo);
                OnsetIndices = find((onsets >= Bouts(BoutNo, 1)) & (onsets <= Bouts(BoutNo, 2)));
                Time = linspace(T(1), T(end), length(m_AM));

                for NoteIndex = 1:length(OnsetIndices),
                    [minval, StartI] = min(abs(Time - onsets(OnsetIndices(NoteIndex))));
                    [minval, EndI] = min(abs(Time - offsets(OnsetIndices(NoteIndex))));
                    Bout.FileName{BoutIndex}{NoteIndex} = SongFile;
                    Bout.AM{BoutIndex}{NoteIndex} = mean(m_AM(StartI:EndI));
                    Bout.FM{BoutIndex}{NoteIndex} = mean(m_FM(StartI:EndI));
                    Bout.Entropy{BoutIndex}{NoteIndex} = mean(m_Entropy(StartI:EndI));
                    Bout.Amp{BoutIndex}{NoteIndex} = mean(m_amplitude(StartI:EndI));
                    Bout.Freq{BoutIndex}{NoteIndex} = mean(m_Freq(StartI:EndI));
                    Bout.PG{BoutIndex}{NoteIndex} = mean(m_PitchGoodness(StartI:EndI));
                    Bout.Pitch{BoutIndex}{NoteIndex} = mean(Pitch_chose(StartI:EndI));
                    Bout.Duration{BoutIndex}{NoteIndex} = Time(EndI) - Time(StartI);
                end
                Bout.SyllDurations{BoutIndex} = offsets(OnsetIndices) - onsets(OnsetIndices);
                Bout.GapDurations{BoutIndex} = onsets(OnsetIndices(2:end)) - offsets(OnsetIndices(1:end-1));
                Bout.BoutLength{BoutIndex} = Bouts(BoutNo,2) - Bouts(BoutNo, 1);
                BoutIndex = BoutIndex + 1;
                fprintf('\t');
            end
        end
        if (strfind(PlotOption, 'on'))
            PlotSpectrogram(DirectoryName, SongFile, FileType);
            hold on;
            plot(SongTime, (Power + (abs(NoiseMSD(1))))*1500);
            for j = 1:length(onsets),
                plot([onsets(j) onsets(j) offsets(j) offsets(j)], [0 8000 8000 0],'r');
            end
            if (~isempty(Bouts))
                for j = 1:size(Bouts,1),
                    plot([Bouts(j,1) Bouts(j,2)], [7500 7500],'b', 'LineWidth', 2);
                end
            end
    %        plot([T(1) T(end)], [((NoiseMSD(1) + ThreshMultiplier*NoiseMSD(2)) + (abs(NoiseMSD(1))))*1500 ((NoiseMSD(1) + ThreshMultiplier*NoiseMSD(2)) + (abs(NoiseMSD(1))))*1500], 'g', 'LineWidth', 2);
            plot([T(1) T(end)], [(Threshold + (abs(NoiseMSD(1))))*1500 (Threshold + (abs(NoiseMSD(1))))*1500], 'c', 'LineWidth', 2);        
            zoom xon;
            uiwait(gcf);
        else
            fprintf('\n');
        end
    catch
        disp(['Could not analyse ', SongFile]);
    end
    clear onsets offsets Bouts BoutOnsets BoutOffsets Temp;
    SongFile = fgetl(Fid);
end
disp(['Discarded ', num2str(DiscardedBouts), ' bouts out of a total of ', num2str(BoutIndex - 1), ' bouts']);
fclose(Fid);
