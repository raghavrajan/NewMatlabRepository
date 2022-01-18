function [FFTLogAmplitudes, BaselineAmpValue, SyllableStatus, Amp_Fs] = LSINA_CalcFFTLogAmpValues(BirdParameters)

% SAPLogAmplitudes - the raw amplitude traces with some padding before and
% after set by the padding variable below
% BaselineAmpValue - for each file, this is obtained by fitting two
% gaussians to the entire log amplitude trace - assumed to be fitting noise
% and sound. The lower mean is the one for the noise and this is taken as
% the baseline amplitude value - BaselineAmpValue
% SyllableStatus - 0 or 1 depending on whether there was enough data for
% padding or not.

Padding = 0.025; % in seconds
WindowSize = 0.005; % in seconds
WindowStep = 0.001; % in seconds
if (BirdParameters.Continuousdata == 0)
    for i = 1:length(BirdParameters.SongFileNames),
        if (isempty(BirdParameters.NoteInfo{i}.onsets))
            FFTLogAmplitudes{i} = [];
            BaselineAmpValue{i} = [];
            SyllableStatus{i} = [];
            Amp_Fs{i} = [];
            continue;
        end
        [RawSong, Fs] = ASSLGetRawData(BirdParameters.DataDirectory, BirdParameters.SongFileNames{i}, BirdParameters.FileType, 0);
        [LogAmplitude] = LSINA_CalcFFTLogAmplitude(RawSong, Fs, WindowSize, WindowStep);
        
        Time = linspace(0, length(RawSong)/Fs, length(LogAmplitude));
        RawSongFs = Fs;
        Fs = 1/(Time(2) - Time(1));
        
        Temp = gmdistribution.fit(LogAmplitude(:), 2);
        
        for j = 1:length(BirdParameters.NoteInfo{i}.onsets),
            SyllableStatus{i}(j) = 1;
            SyllableOnsetTimeIndex = round(((BirdParameters.NoteInfo{i}.onsets(j)/1000) - Padding) * Fs);
            SyllableOffsetTimeIndex = round(((BirdParameters.NoteInfo{i}.offsets(j)/1000) + Padding) * Fs);
            if (SyllableOnsetTimeIndex(1) <= 0)
                SyllableOnsetTimeIndex(1) = 1;
                SyllableStatus{i}(j) = 0;
            end
            if (SyllableOffsetTimeIndex > length(LogAmplitude))
                SyllableOffsetTimeIndex = length(LogAmplitude);
                SyllableStatus{i}(j) = 0;
            end
            FFTLogAmplitudes{i}{j} = (LogAmplitude(SyllableOnsetTimeIndex:SyllableOffsetTimeIndex));
            BaselineAmpValue{i}(j) = min(Temp.mu);
            Amp_Fs{i}(j) = Fs;
        end
    end
end
