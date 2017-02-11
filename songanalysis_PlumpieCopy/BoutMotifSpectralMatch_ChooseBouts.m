function [BoutLimits] = BoutMotifSpectralMatch_ChooseBouts(DirectoryName, FileList, FileType, ChooseBouts, GetRidNoise)

InterBoutInterval = 1;
MinBoutLen = 1;

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

Fid = fopen(FileList, 'r');
Temp = textscan(Fid, '%s', 'delimiter', '\n');
fclose(Fid);

SongFiles = Temp{1};

cd(DirectoryName);
for SongFileNo = 1:length(SongFiles),
    SongFile = SongFiles{SongFileNo};
    if (isempty(SongFile))
        continue;
    end
    
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

   try
    
    SongTime = (1:1:length(Song))/Fs;
    
    if (strfind(ChooseBouts, 'FullFile'))
        Bouts = [SongTime(1) SongTime(end)];
        BoutLimits{SongFileNo} = Bouts;
        clear Bouts;
        continue;
    end
    
    WinSize = round(FFTWinSize * Fs);
    WinOverlap = round(FFTWinOverlap * WinSize);

    [S1, F, T, P] = spectrogram(Song, hamming(WinSize), WinOverlap, WinSize, Fs);
    
    Freq1 = find((F >= 860) & (F <= 8600));
    Power = log10(sum(S1(Freq1,:).*conj(S1(Freq1,:))));
    S = log10(abs(S1(Freq1,:)));
 
    clear Obj;
    disp(['Removed ', num2str(length(find(isinf(Power)))), ' values from the variable Power (length = ', num2str(length(Power)), ')']);
    
    Power(isinf(Power)) = min(Power(~isinf(Power)));
    
    Obj = gmdistribution.fit(Power', 2, 'Start', 'randSample', 'Replicates', 10);
    
    [SoundMSD(1), Index] = max(Obj.mu);
    SoundMSD(2) = Obj.Sigma(Index);
    [NoiseMSD(1), Index] = min(Obj.mu);
    NoiseMSD(2) = Obj.Sigma(Index);

    Threshold = mean(Obj.mu);
    disp(['Threshold is ', num2str(Threshold)]);
    
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
        if (strfind(ChooseBouts, 'yes'))
            PlotSpectrogram(DirectoryName, SongFile, FileType);
            zoom xon;
            Flag = 1;
            Bouts = [];
            while(Flag)
                FileYesOrNo = inputdlg('Do you want to use this file? (y for yes and n for no)');
                if (isempty(strfind(FileYesOrNo{1}, 'y')))
                    FileYesNo = 0;
                    break;
                else
                    FileYesNo = 1;
                end
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
            if (FileYesNo == 0)
                clear onsets offsets Bouts BoutOnsets BoutOffsets Temp;
                continue;
            end
        else
            if (strfind(ChooseBouts, 'FullFile'))
                Bouts = [SongTime(1) SongTime(end)];
            end
        end
    end
    BoutLimits{SongFileNo} = Bouts;
    clear Bouts;
   end
end
disp('Finished')