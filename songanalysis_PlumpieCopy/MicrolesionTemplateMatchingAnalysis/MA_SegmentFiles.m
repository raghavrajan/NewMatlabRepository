function [onsets, offsets, durs, gaps, threshold, len, Fs] = MA_SegmentFiles(SongFileName, FileType, DataDir, OutputDir)

% Common variables
min_int = 7; % minimum interval between syllables in ms
min_dur = 7; % minimum syllable duration in ms

OutputFile = [OutputDir, SongFileName, '.not.mat'];

[Song, Fs] = MA_ReadSongFile(DataDir, SongFileName, FileType);

if (exist(OutputFile, 'file'))
    load(OutputFile);
else
    [LogAmplitude] = MA_CalculateLogAmplitudeAronovFee(Song, Fs);
    
    threshold = ASSLCalculateFisherThreshold(LogAmplitude);
    [onsets, offsets] = ASSLSegmentDataAronovFee(LogAmplitude, Fs, min_int, min_dur, threshold);
    labels = char(ones(1, length(onsets))*double('0'));
    sm_win = 2.5;
    save(OutputFile, 'onsets', 'offsets', 'labels', 'min_int', 'min_dur', 'sm_win', 'threshold');
end

% Calculated durations (durs) and gaps and the duration of the entire file
% (len)
durs = offsets - onsets;
gaps = onsets(2:end) - offsets(1:end-1);
len = length(Song) * 1000/Fs;