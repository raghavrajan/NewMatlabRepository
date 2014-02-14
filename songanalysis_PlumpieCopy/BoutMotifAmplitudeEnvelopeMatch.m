function [Bout] = BoutMotifAmplitudeEnvelopeMatch(DirectoryName, FileList, FileType, PlotOption, ThreshMultiplier, MotifTemplate, ChooseBouts, StretchValues, GetRidNoise)

InterBoutInterval = 1;
MinBoutLen = 3;

TimeNow = datestr(now, 'mmddyyHHMMSS');

% Open log file for results
if (ispc)
    LogFid = fopen(['C:\Documents and Settings\PlumPie\My Documents\RaghavData\BoutStats_',TimeNow,'.log'], 'w');
else
    LogFid = fopen(['/home/raghav/BoutStatistics/BoutStats_',TimeNow,'.log'], 'w');
end

% Write input parameters to Log file
fprintf(LogFid, 'DirectoryName: %s\n', DirectoryName);
fprintf(LogFid, 'FileList: %s\n', FileList);
fprintf(LogFid, 'FileType: %s\n', FileType);
fprintf(LogFid, 'PlotOption: %s\n', PlotOption);
fprintf(LogFid, 'FileType: %s\n', FileType);
fprintf(LogFid, 'Threshold multiplier is %i\n', ThreshMultiplier);
fprintf(LogFid, 'Choose bouts: %s\n', ChooseBouts);
fprintf(LogFid, 'Stretch values are ');
for i = 1:length(StretchValues),
    fprintf(LogFid, '%i,', StretchValues(i));
end
fprintf(LogFid, '\n\n\n');

min_int = 0.005;
min_dur = 0.010;
sm_win = 0.008;

FFTWinSize = sm_win; % in sec
FFTWinOverlap = 0.5;

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

x = 1:1:size(MotifTemplate, 2);

BoutIndex = 1;

Fid = fopen(FileList, 'r');
SongFile = fgetl(Fid);

cd(DirectoryName);
while (ischar(SongFile(1)))
    Slash = find((SongFile == '/') | (SongFile == '\'));
    if (~isempty(Slash))
        SongFile = SongFile(Slash(end)+1:end);
    end
    disp(SongFile);
    fprintf(LogFid, '\nSongfile name: %s\n', SongFile);
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
   % [m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude , m_Freq, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight ]=deriv(Song, Fs);

   try
    SongTime = (1:1:length(Song))/Fs;
    WinSize = round(FFTWinSize * Fs);
    WinOverlap = round(FFTWinOverlap * WinSize);

    [S1, F, T, P] = spectrogram(Song, hamming(WinSize), WinOverlap, WinSize, Fs);
    
    Freq1 = find((F >= 860) & (F <= 8600));
    Power = log10(sum(S1(Freq1,:).*conj(S1(Freq1,:))));
    S = log10(abs(S1(Freq1,:)));
    
    Obj = gmdistribution.fit(Power', 2, 'Start', 'randSample', 'Replicates', 10);
    [SoundMSD(1), Index] = max(Obj.mu);
    SoundMSD(2) = Obj.Sigma(Index);
    [NoiseMSD(1), Index] = min(Obj.mu);
    NoiseMSD(2) = Obj.Sigma(Index);

    Threshold = mean(Obj.mu);
    disp(['Threshold is ', num2str(Threshold)]);
    
    fprintf(LogFid, 'Noise mean and var are %g and %g\n', NoiseMSD(1), NoiseMSD(2));
    fprintf(LogFid, 'Sound mean and var are %g and %g\n', SoundMSD(1), SoundMSD(2));
    fprintf(LogFid, 'Threshold is %g\n', Threshold);
    
    % Power = 10*log10(sum(P));
    % ThreshCrossing = Power > Threshold;
    Power = spline(T, Power, SongTime);
    [onsets, offsets] = segment_song_for_labelling(Power, Fs, min_int*1000, min_dur*1000, Threshold, NoiseMSD);
    onsets = onsets/1000;
    offsets = offsets/1000;

    if (size(onsets,1) < size(onsets, 2))
        onsets = onsets';
    end
    if (size(offsets,1) < size(offsets, 2))
        offsets = offsets';
    end
    
    if (strfind(GetRidNoise, 'yes'))
            NoiseFreq = find(F < 500);
            NoisePower = log10(sum(S1(NoiseFreq,:).*conj(S1(NoiseFreq,:))));
            Width = 0.02;
            XGauss = 1:1:(1 + round(2 * 4 * Width * (1/(T(2) - T(1)))));
            XGauss = XGauss - (length(XGauss) + 1)/2;
            GaussWin = (1/((Width * (1/(T(2) - T(1)))) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * (1/(T(2) - T(1)))) * (Width * (1/(T(2) - T(1))))));
            NTC = NoisePower > 1.9;
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
    
    if (strfind(ChooseBouts, 'no'))

        Temp = (onsets(2:end) - offsets(1:end-1));
        if (onsets(1) > InterBoutInterval)
            BoutOnsets = [onsets(1); onsets(find(Temp > InterBoutInterval) + 1)];
            if (offsets(end) > (T(end) - InterBoutInterval))
                BoutOffsets = [offsets(Temp > InterBoutInterval); (offsets(end) - InterBoutInterval)];
            else
                BoutOffsets = [offsets(Temp > InterBoutInterval); offsets(end)];
            end
        else
            BoutOnsets = [(onsets(1) + InterBoutInterval); onsets(find(Temp > InterBoutInterval) + 1)];
            if (offsets(end) > (T(end) - InterBoutInterval))
                BoutOffsets = [offsets(Temp > InterBoutInterval); (offsets(end) - InterBoutInterval)];
            else
                BoutOffsets = [offsets(Temp > InterBoutInterval); offsets(end)];
            end
        end
        BoutOnsets = BoutOnsets - InterBoutInterval;
        BoutOffsets = BoutOffsets + InterBoutInterval;

        Bouts = [BoutOnsets BoutOffsets];
    else
        if (strfind('ChooseBouts', 'yes'))
            PlotSpectrogram(DirectoryName, SongFile, FileType);
            zoom xon;
            Flag = 1;
            Bouts = [];
            while(Flag)
                TempMsgBox = msgbox('Zoom into the desired region and then close this box when ready');
                waitfor(TempMsgBox);
                disp('Choose bout onset and offset');
                [x1] = ginput(2);
                Bouts(end+1,:) = [x1(1) x1(2)];
                YesNo = inputdlg('Do you want to choose more bouts? (y for yes and n for no)');
                if (isempty(strfind(YesNo{1}, 'y')))
                    Flag = 0;
                end
            end
        else
            if (strfind(ChooseBouts, 'FullFile'))
                Bouts = [SongTime(1) SongTime(end)];
            end
        end
    end
    if (~isempty(Bouts))
        DiscardedBouts = DiscardedBouts + length(find((Bouts(:,2) - Bouts(:,1)) <= MinBoutLen));
        Bouts = Bouts((Bouts(:,2) - Bouts(:,1)) > MinBoutLen,:);

%          BoutRemoval = [];
%          for BoutNo = 1:size(Bouts,1),
%              OnsetIndices = find((onsets >= Bouts(BoutNo, 1)) & (onsets <= Bouts(BoutNo, 2)));
%              ISIs = onsets(OnsetIndices(2:end)) - offsets(OnsetIndices(1:end-1));
%              if ((length(find(ISIs > 0.4))/length(ISIs)) >= 0.1)
%                  BoutRemoval = [BoutRemoval; BoutNo];
%              end
%          end
%          DiscardedBouts = DiscardedBouts + length(BoutRemoval);
%          if (length(BoutRemoval) > 0)
%              Bouts(BoutRemoval,:) = [];
%          end
        %Time = linspace(T(1), T(end), length(m_Entropy));
        for BoutNo = 1:size(Bouts,1),
            fprintf(LogFid, 'Bout #%i\n', BoutNo);
            fprintf('B%i ', BoutNo);
            Bout.BoutAmplitudes{BoutIndex} = Power((T >= Bouts(BoutNo,1)) & (T <= Bouts(BoutNo, 2)));
         %   Bout.BoutPG{BoutIndex} = m_PitchGoodness(find((Time >= Bouts(BoutNo,1)) & (Time <= Bouts(BoutNo, 2))));
          %  Bout.BoutFM{BoutIndex} = m_FM(find((Time >= Bouts(BoutNo,1)) & (Time <= Bouts(BoutNo, 2))));
          %  Bout.BoutEntropy{BoutIndex} = m_Entropy(find((Time >= Bouts(BoutNo,1)) & (Time <= Bouts(BoutNo, 2))));
            Bout.BoutAmplitudeFs{BoutIndex} = 1/(T(2) - T(1));
         %   Bout.BoutPGFs{BoutIndex} = 1/(Time(2) - Time(1));
            Indices = (find((T >= Bouts(BoutNo,1)) & (T <= Bouts(BoutNo, 2))));

            WarpIndex = 1;
            % Matches = [];
            for i = [StretchValues],
                xx = linspace(x(1), x(end), (length(MotifTemplate) * (1 + i/100)));
                clear WMotif;
                WMotif = spline(x, MotifTemplate, xx);
                TempS = sum(S(:,Indices));

%                 WinMean = zeros((length(TempS) - length(WMotif) + 1), 1);
%                 WinSTD = zeros((length(TempS) - length(WMotif) + 1), 1);

                for ColNo = 1:(length(TempS) - length(WMotif) + 1),
                    StartIndex = ColNo;
                    WinIndices = StartIndex:1:(StartIndex + length(WMotif) - 1);
%                     WinMean(ColNo) = mean(TempS(WinIndices));
%                     WinSTD(ColNo) = std(TempS(WinIndices));
                    Data = TempS(WinIndices);
                    Data = (Data - mean(Data))/std(Data);
                    Match(ColNo) = 1/sum(sum((abs(Data - WMotif))));
                end
                %[Match] = CalTemplateMatch(WMotif, TempS, WinMean, WinSTD);
                fprintf(LogFid, '%g\n', max(Match));
                Bout.BoutSeqMatch{BoutIndex}{WarpIndex} = Match;
                WarpIndex = WarpIndex + 1;
                fprintf('>');
            end
            Bout.T{BoutIndex} = T(Indices);
            clear Match;
            for MatchNo = 1:length(Bout.BoutSeqMatch{BoutIndex}),
                Match(MatchNo,:) = Bout.BoutSeqMatch{BoutIndex}{MatchNo}(1:length(Bout.BoutSeqMatch{BoutIndex}{end}));
            end
            Match = max(Match);
            Bout.MaxBoutSeqMatch{BoutIndex} = Match;
            clear Match;
            [MaxVal, MaxInd] = max(Bout.MaxBoutSeqMatch{BoutIndex});
            Bout.MaxBoutSeqMatchVal{BoutIndex} = [Bout.T{BoutIndex}(MaxInd) MaxVal];
            fprintf(LogFid, 'Maximum match value is %g at time %g\n', MaxVal, Bout.T{BoutIndex}(MaxInd));
            BoutIndex = BoutIndex + 1;
            fprintf('\t');
        end
    end
    
    if (strfind(PlotOption, 'on'))
        PlotSpectrogram(DirectoryName, SongFile, FileType);
        hold on;
        plot(SongTime, (Power + (abs(min(Power))))*100);
        for j = 1:length(onsets),
            plot([onsets(j) onsets(j) offsets(j) offsets(j)], [0 8000 8000 0],'r');
        end
        if (~isempty(Bouts))
            for j = 1:size(Bouts,1),
                plot([Bouts(j,1) Bouts(j,2)], [7500 7500],'b', 'LineWidth', 2);
                %plot(Bout.T{BoutIndex - j}(1:length(Bout.BoutSeqMatch{BoutIndex - j})), (Bout.BoutSeqMatch{BoutIndex-j}/max(Bout.BoutSeqMatch{BoutIndex-j})) *8000, 'r', 'LineWidth', 2);
            end
        end
        plot([T(1) T(end)], [((NoiseMSD(1) + ThreshMultiplier*NoiseMSD(2)) + (abs(min(Power))))*100 ((NoiseMSD(1) + ThreshMultiplier*NoiseMSD(2)) + (abs(min(Power))))*100], 'g', 'LineWidth', 2);
        plot([T(1) T(end)], [(Threshold + (abs(min(Power))))*100 (Threshold + (abs(min(Power))))*100], 'c', 'LineWidth', 2);        
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
fclose(LogFid);
