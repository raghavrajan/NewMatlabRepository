function [BoutFFT] = MA_AnalyseSavedData_RhythmSpectrogram(Parameters, TitleString, varargin)

if (nargin > 2)
    BirdIndices = varargin{1};
else
    BirdIndices = 1:1:length(Parameters);
end

OutputDir = '/home/raghav/HVC_MicrolesionDataFigures/PaperFigures/';

PrePostDays = [1 2; 2 3; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 2 5; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];

Parameters = Parameters(BirdIndices);
PrePostDays = PrePostDays(BirdIndices, :);

Params.DiscreteBout_Fs = 500; % sampling rate for representation of bouts as 1s and 0s - Hz
Params.NoofBins = 50;
Params.DurEdges = linspace(-2.5, 0, Params.NoofBins);
Params.NFFT = 2^nextpow2(round(2*Params.DiscreteBout_Fs));
Params.Freq = Params.DiscreteBout_Fs/2*linspace(0, 1, Params.NFFT/2+1);
Params.FreqIndices = find(Params.Freq <= 30);
% Params.HanningWinLen = 0.01; % in sec
% Params.HanningWin = hanning(round(Params.HanningWinLen * Parameters.PreDirFs{1}{1}));

for ParameterNo = 1:length(Parameters),
    fprintf('%s >> ', Parameters(ParameterNo).BirdName);
    for i = 1:Parameters(ParameterNo).NoPreDays,
        fprintf('Pre Day #%i >> ', i);
        fprintf('Dir >> ');
        BoutIndex = 1;
        
        for j = 1:length(Parameters(ParameterNo).PreDirBoutLens{i}),
            [RawData, Fs] = ASSLGetRawData(Parameters(ParameterNo).PreDataDir{i}, Parameters(ParameterNo).PreDirSongFileNames{i}{j}, Parameters(ParameterNo).FileType, 0);
            for k = 1:length(Parameters(ParameterNo).PreDirBoutLens{i}{j}),
                SongBout = RawData(ceil((Parameters(ParameterNo).PreDirBoutOnsets{i}{j}(k) - 0) * Fs/1000):floor((Parameters(ParameterNo).PreDirBoutOffsets{i}{j}(k) + 0)*Fs/1000));
                
                [TempBoutFFT] = MA_CalcRhythmSpectrogram(SongBout, Fs, Params.DiscreteBout_Fs, Params.Freq(Params.FreqIndices));
                BoutFFT(ParameterNo).Dir{i}(BoutIndex,:) = TempBoutFFT;
                
                BoutIndex = BoutIndex + 1;
            end
        end

        fprintf('Undir >> ');
        
        BoutIndex = 1;
        
        for j = 1:length(Parameters(ParameterNo).PreUnDirBoutLens{i}),
            [RawData, Fs] = ASSLGetRawData(Parameters(ParameterNo).PreDataDir{i}, Parameters(ParameterNo).PreUnDirSongFileNames{i}{j}, Parameters(ParameterNo).FileType, 0);
            for k = 1:length(Parameters(ParameterNo).PreUnDirBoutLens{i}{j}),
                SongBout = RawData(ceil((Parameters(ParameterNo).PreUnDirBoutOnsets{i}{j}(k) - 0) * Fs/1000):floor((Parameters(ParameterNo).PreUnDirBoutOffsets{i}{j}(k) + 0)*Fs/1000));
                
                [TempBoutFFT] = MA_CalcRhythmSpectrogram(SongBout, Fs, Params.DiscreteBout_Fs, Params.Freq(Params.FreqIndices));
                BoutFFT(ParameterNo).UnDir{i}(BoutIndex,:) = TempBoutFFT;
                BoutIndex = BoutIndex + 1;
            end
        end
    end


    for i = 1:Parameters(ParameterNo).NoPostDays,
        fprintf('Post Day #%i >> ', i);
        fprintf('Dir >> ');
        BoutIndex = 1;
        for j = 1:length(Parameters(ParameterNo).PostDirBoutLens{i}),
            [RawData, Fs] = ASSLGetRawData(Parameters(ParameterNo).PostDataDir{i}, Parameters(ParameterNo).PostDirSongFileNames{i}{j}, Parameters(ParameterNo).FileType, 0);
            for k = 1:length(Parameters(ParameterNo).PostDirBoutLens{i}{j}),
                SongBout = RawData(ceil((Parameters(ParameterNo).PostDirBoutOnsets{i}{j}(k) - 0) * Fs/1000):floor((Parameters(ParameterNo).PostDirBoutOffsets{i}{j}(k) + 0)*Fs/1000));
                
                [TempBoutFFT] = MA_CalcRhythmSpectrogram(SongBout, Fs, Params.DiscreteBout_Fs, Params.Freq(Params.FreqIndices));
                BoutFFT(ParameterNo).Dir{i + Parameters(ParameterNo).NoPreDays}(BoutIndex,:) = TempBoutFFT;
                
                BoutIndex = BoutIndex + 1;
            end
        end

        fprintf('Undir >> ');
        BoutIndex = 1;
        for j = 1:length(Parameters(ParameterNo).PostUnDirBoutLens{i}),
            [RawData, Fs] = ASSLGetRawData(Parameters(ParameterNo).PostDataDir{i}, Parameters(ParameterNo).PostUnDirSongFileNames{i}{j}, Parameters(ParameterNo).FileType, 0);
            for k = 1:length(Parameters(ParameterNo).PostUnDirBoutLens{i}{j}),
                SongBout = RawData(ceil((Parameters(ParameterNo).PostUnDirBoutOnsets{i}{j}(k) - 0) * Fs/1000):floor((Parameters(ParameterNo).PostUnDirBoutOffsets{i}{j}(k) + 0)*Fs/1000));
                
                [TempBoutFFT] = MA_CalcRhythmSpectrogram(SongBout, Fs, Params.DiscreteBout_Fs, Params.Freq(Params.FreqIndices));
                BoutFFT(ParameterNo).UnDir{i + Parameters(ParameterNo).NoPreDays}(BoutIndex,:) = TempBoutFFT;
                
                BoutIndex = BoutIndex + 1;
            end
        end

    end
    fprintf('\n');
end

Colors = 'rgbcmk';
for i = 1:length(BoutFFT),
    figure;
    for j = 1:2,
        plot(Params.Freq(Params.FreqIndices), mean(BoutFFT(i).Dir{PrePostDays(i,j)}), Colors(j), 'LineWidth', 2);
        hold on;
    end
    
    for j = 1:2,
        plot(Params.Freq(Params.FreqIndices), mean(BoutFFT(i).UnDir{PrePostDays(i,j)}), Colors(j+2), 'LineWidth', 2);
        hold on;
    end
    title(Parameters(i).BirdName, 'FontSize', 16, 'FontWeight', 'bold');
end
        
% Now to calculate different parameters related to the rhythm spectrogram
% based on Goldberg and Fee (2011) and Kaie and Hahnloser (2011)

Index = find(Params.Freq(Params.FreqIndices) >= 3, 1, 'first');

%==========================================================================
% Plot the value of the max of the peak between 3 and 30Hz

clear DirMedians UnDirMedians;
for i = 1:length(BoutFFT),
    TempBoutFFT = mean(BoutFFT(i).Dir{PrePostDays(i,1)});
    DirMedians(i,1) = max(TempBoutFFT(Index:end));
    
    TempBoutFFT = mean(BoutFFT(i).Dir{PrePostDays(i,2)});
    DirMedians(i,2) = max(TempBoutFFT(Index:end));
    
    TempBoutFFT = mean(BoutFFT(i).UnDir{PrePostDays(i,1)});
    UnDirMedians(i,1) = max(TempBoutFFT(Index:end));
    
    TempBoutFFT = mean(BoutFFT(i).UnDir{PrePostDays(i,2)});
    UnDirMedians(i,2) = max(TempBoutFFT(Index:end));
end

MA_PlotVsLesionSize([Parameters.PercentTotalHVCremaining], DirMedians, UnDirMedians, ['Peak of FFT (Post/Pre)'], OutputDir, 'NorthWest', ['FFTPeak'], TitleString);

MA_PlotPreVsPost(DirMedians, UnDirMedians, ['Peak of FFT'], OutputDir, ['FFTPeak'], TitleString);

MA_PlotDirVsUnDir([DirMedians(:,1) UnDirMedians(:,1)], [DirMedians(:,2) UnDirMedians(:,2)], ['Peak of FFT'], OutputDir, ['FFTPeak'], TitleString);

%==========================================================================
% Plot the modulation between 3 and 30Hz

    LinearFit_Coeffs = polyfit(Params.Freq(Params.FreqIndices(Index:end)), TempBoutFFT(Index:end), 1);
    LinearFit = polyval(LinearFit_Coeffs, Params.Freq(Params.FreqIndices(Index:end)));
    Results.DirRhythm_Modulation(i) = sqrt(sum((LinearFit - TempBoutFFT(Index:end)).^2));

clear DirMedians UnDirMedians;
for i = 1:length(BoutFFT),
    TempBoutFFT = mean(BoutFFT(i).Dir{PrePostDays(i,1)});
    LinearFit_Coeffs = polyfit(Params.Freq(Params.FreqIndices(Index:end)), TempBoutFFT(Index:end), 1);
    LinearFit = polyval(LinearFit_Coeffs, Params.Freq(Params.FreqIndices(Index:end)));
    DirMedians(i,1) = sqrt(sum((LinearFit - TempBoutFFT(Index:end)).^2));
    
    TempBoutFFT = mean(BoutFFT(i).Dir{PrePostDays(i,2)});
    LinearFit_Coeffs = polyfit(Params.Freq(Params.FreqIndices(Index:end)), TempBoutFFT(Index:end), 1);
    LinearFit = polyval(LinearFit_Coeffs, Params.Freq(Params.FreqIndices(Index:end)));
    DirMedians(i,2) = sqrt(sum((LinearFit - TempBoutFFT(Index:end)).^2));
    
    TempBoutFFT = mean(BoutFFT(i).UnDir{PrePostDays(i,1)});
    LinearFit_Coeffs = polyfit(Params.Freq(Params.FreqIndices(Index:end)), TempBoutFFT(Index:end), 1);
    LinearFit = polyval(LinearFit_Coeffs, Params.Freq(Params.FreqIndices(Index:end)));
    UnDirMedians(i,1) = sqrt(sum((LinearFit - TempBoutFFT(Index:end)).^2));
    
    TempBoutFFT = mean(BoutFFT(i).UnDir{PrePostDays(i,2)});
    LinearFit_Coeffs = polyfit(Params.Freq(Params.FreqIndices(Index:end)), TempBoutFFT(Index:end), 1);
    LinearFit = polyval(LinearFit_Coeffs, Params.Freq(Params.FreqIndices(Index:end)));
    UnDirMedians(i,2) = sqrt(sum((LinearFit - TempBoutFFT(Index:end)).^2));
end

MA_PlotVsLesionSize([Parameters.PercentTotalHVCremaining], DirMedians, UnDirMedians, ['Modulation of FFT (Post/Pre)'], OutputDir, 'NorthWest', ['FFTModulation'], TitleString);

MA_PlotPreVsPost(DirMedians, UnDirMedians, ['Modulation of FFT'], OutputDir, ['FFTModulation'], TitleString);

MA_PlotDirVsUnDir([DirMedians(:,1) UnDirMedians(:,1)], [DirMedians(:,2) UnDirMedians(:,2)], ['Modulation of FFT'], OutputDir, ['FFTModulation'], TitleString);

disp('Finished analysing rhythm spectrogram');