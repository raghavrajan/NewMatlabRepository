function [] = MA_AnalyseSavedData_SpectralCorr(Parameters, TitleString, varargin)

if (nargin > 2)
    BirdIndices = varargin{1};
else
    BirdIndices = 1:1:length(Parameters);
end

OutputDir = '/home/raghav/HVC_MicrolesionDataFigures/PaperFigures/';

PrePostDays = [1 2; 2 3; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2; 2 5; 1 2; 1 2; 1 2; 1 2; 1 2; 1 2];

Parameters = Parameters(BirdIndices);
PrePostDays = PrePostDays(BirdIndices, :);

WindowSize = 10; % window size in ms for calculating PSD
StepSize = 1; % step in ms for calculating PSD
BandWidth = 1.5; % bandwith for generating multi-taper spectrogram
BoutSize = 3; % minimum length of bout to be considered in sec

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
                
                LogAmplitude = ASSLCalculateLogAmplitudeAronovFee(SongBout, Fs);
                BoutAmplitude(ParameterNo).Dir{i}{BoutIndex} = spline((1:1:length(LogAmplitude))/Fs, LogAmplitude, (1:1:length(LogAmplitude)*Params.DiscreteBout_Fs/Fs)/Params.DiscreteBout_Fs);
                
                BoutIndex = BoutIndex + 1;
            end
        end

        fprintf('Undir >> ');
        
        BoutIndex = 1;
        
        for j = 1:length(Parameters(ParameterNo).PreUnDirBoutLens{i}),
            [RawData, Fs] = ASSLGetRawData(Parameters(ParameterNo).PreDataDir{i}, Parameters(ParameterNo).PreUnDirSongFileNames{i}{j}, Parameters(ParameterNo).FileType, 0);
            for k = 1:length(Parameters(ParameterNo).PreUnDirBoutLens{i}{j}),
                SongBout = RawData(ceil((Parameters(ParameterNo).PreUnDirBoutOnsets{i}{j}(k) - 0) * Fs/1000):floor((Parameters(ParameterNo).PreUnDirBoutOffsets{i}{j}(k) + 0)*Fs/1000));
                
                LogAmplitude = ASSLCalculateLogAmplitudeAronovFee(SongBout, Fs);
                BoutAmplitude(ParameterNo).UnDir{i}{BoutIndex} = spline((1:1:length(LogAmplitude))/Fs, LogAmplitude, (1:1:length(LogAmplitude)*Params.DiscreteBout_Fs/Fs)/Params.DiscreteBout_Fs);
                
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
                
                LogAmplitude = ASSLCalculateLogAmplitudeAronovFee(SongBout, Fs);
                BoutAmplitude(ParameterNo).Dir{i + Parameters(ParameterNo).NoPreDays}{BoutIndex} = spline((1:1:length(LogAmplitude))/Fs, LogAmplitude, (1:1:length(LogAmplitude)*Params.DiscreteBout_Fs/Fs)/Params.DiscreteBout_Fs);
                
                BoutIndex = BoutIndex + 1;
            end
        end

        fprintf('Undir >> ');
        BoutIndex = 1;
        for j = 1:length(Parameters(ParameterNo).PostUnDirBoutLens{i}),
            [RawData, Fs] = ASSLGetRawData(Parameters(ParameterNo).PostDataDir{i}, Parameters(ParameterNo).PostUnDirSongFileNames{i}{j}, Parameters(ParameterNo).FileType, 0);
            for k = 1:length(Parameters(ParameterNo).PostUnDirBoutLens{i}{j}),
                SongBout = RawData(ceil((Parameters(ParameterNo).PostUnDirBoutOnsets{i}{j}(k) - 0) * Fs/1000):floor((Parameters(ParameterNo).PostUnDirBoutOffsets{i}{j}(k) + 0)*Fs/1000));
                
                LogAmplitude = ASSLCalculateLogAmplitudeAronovFee(SongBout, Fs);
                BoutAmplitude(ParameterNo).UnDir{i + Parameters(ParameterNo).NoPreDays}{BoutIndex} = spline((1:1:length(LogAmplitude))/Fs, LogAmplitude, (1:1:length(LogAmplitude)*Params.DiscreteBout_Fs/Fs)/Params.DiscreteBout_Fs);
                
                BoutIndex = BoutIndex + 1;
            end
        end

    end
    fprintf('\n');
end


%============= Now do the spectral correlations ===========================
% First pick 10 random bouts for each day and each condition with bout
% length > 3 seconds.

WindowSize = 10; % window size in ms for calculating PSD
StepSize = 1; % step in ms for calculating PSD
BandWidth = 1.5; % bandwith for generating multi-taper spectrogram
BoutSize = 3; % minimum length of bout to be considered in sec

disp('Doing spectral correlations ...');
% First for pre song

% First for undir song

FileTypeCellArray = cellstr(char(ones(length(Parameters.PreUnDirBouts), 1) * double(Parameters.FileType)))';
WindowSizeCellArray = num2cell(ones(1, length(Parameters.PreUnDirBouts)) * WindowSize);
StepSizeCellArray = num2cell(ones(1, length(Parameters.PreUnDirBouts)) * StepSize);
BandWidthCellArray = num2cell(ones(1, length(Parameters.PreUnDirBouts)) * BandWidth);
BoutSizeCellArray = num2cell(ones(1, length(Parameters.PreUnDirBouts)) * BoutSize);

[Parameters.PreUnDirSpectralCorrs] = cellfun(@MA_DoSpectralCorrelations, Parameters.PreUnDirBouts, Parameters.PreUnDirSongFileNames, Parameters.PreDataDir, FileTypeCellArray, WindowSizeCellArray, StepSizeCellArray, BandWidthCellArray, BoutSizeCellArray, 'UniformOutput', 0);

% Now for post song
% First for undir song

FileTypeCellArray = cellstr(char(ones(length(Parameters.PostUnDirBouts), 1) * double(Parameters.FileType)))';
WindowSizeCellArray = num2cell(ones(1, length(Parameters.PostUnDirBouts)) * WindowSize);
StepSizeCellArray = num2cell(ones(1, length(Parameters.PostUnDirBouts)) * StepSize);
BandWidthCellArray = num2cell(ones(1, length(Parameters.PostUnDirBouts)) * BandWidth);
BoutSizeCellArray = num2cell(ones(1, length(Parameters.PostUnDirBouts)) * BoutSize);

[Parameters.PostUnDirSpectralCorrs] = cellfun(@MA_DoSpectralCorrelations, Parameters.PostUnDirBouts, Parameters.PostUnDirSongFileNames, Parameters.PostDataDir, FileTypeCellArray, WindowSizeCellArray, StepSizeCellArray, BandWidthCellArray, BoutSizeCellArray, 'UniformOutput', 0);

% First for dir song

FileTypeCellArray = cellstr(char(ones(length(Parameters.PreDirBouts), 1) * double(Parameters.FileType)))';
WindowSizeCellArray = num2cell(ones(1, length(Parameters.PreDirBouts)) * WindowSize);
StepSizeCellArray = num2cell(ones(1, length(Parameters.PreDirBouts)) * StepSize);
BandWidthCellArray = num2cell(ones(1, length(Parameters.PreDirBouts)) * BandWidth);
BoutSizeCellArray = num2cell(ones(1, length(Parameters.PreDirBouts)) * BoutSize);

[Parameters.PreDirSpectralCorrs] = cellfun(@MA_DoSpectralCorrelations, Parameters.PreDirBouts, Parameters.PreDirSongFileNames, Parameters.PreDataDir, FileTypeCellArray, WindowSizeCellArray, StepSizeCellArray, BandWidthCellArray, BoutSizeCellArray, 'UniformOutput', 0);

% Now for post song
% First for undir song

FileTypeCellArray = cellstr(char(ones(length(Parameters.PostDirBouts), 1) * double(Parameters.FileType)))';
WindowSizeCellArray = num2cell(ones(1, length(Parameters.PostDirBouts)) * WindowSize);
StepSizeCellArray = num2cell(ones(1, length(Parameters.PostDirBouts)) * StepSize);
BandWidthCellArray = num2cell(ones(1, length(Parameters.PostDirBouts)) * BandWidth);
BoutSizeCellArray = num2cell(ones(1, length(Parameters.PostDirBouts)) * BoutSize);

[Parameters.PostDirSpectralCorrs] = cellfun(@MA_DoSpectralCorrelations, Parameters.PostDirBouts, Parameters.PostDirSongFileNames, Parameters.PostDataDir, FileTypeCellArray, WindowSizeCellArray, StepSizeCellArray, BandWidthCellArray, BoutSizeCellArray, 'UniformOutput', 0);
%==========================================================================