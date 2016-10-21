function [onsets, offsets, durs, gaps, threshold, len, Fs] = MA_SegmentFiles(SongFileName, FileType, DataDir, OutputDir, BoutOnsets, BoutOffsets, BoutOrNot)

% Common variables
min_int = 7; % minimum interval between syllables in ms
min_dur = 7; % minimum syllable duration in ms

OutputFile = [OutputDir, SongFileName, '.not.mat'];

[Song, Fs] = MA_ReadSongFile(DataDir, SongFileName, FileType);

if (exist(OutputFile, 'file'))
    load(OutputFile);
else
    labels = [];
    onsets = [];
    offsets = [];
    if (~isempty(BoutOnsets))
        for i = 1:length(BoutOnsets),
            [LogAmplitude] = MA_CalculateLogAmplitudeAronovFee(Song(round(BoutOnsets(i)*Fs/1000):round(BoutOffsets(i)*Fs/1000)), Fs);

            threshold(i,:) = ASSLCalculateFisherThreshold(LogAmplitude);
            [temp_onsets, temp_offsets] = ASSLSegmentDataAronovFee(LogAmplitude, Fs, min_int, min_dur, threshold(i,:), BoutOrNot);
            labels = [labels char(ones(1, length(temp_onsets))*double('0'))];
            if (~isempty(temp_onsets))
                onsets = [onsets; (temp_onsets(:) + BoutOnsets(i))];
                offsets = [offsets; (temp_offsets(:) + BoutOnsets(i))];
            end
        end    
    else
        threshold = [];
    end
    sm_win = 2.5;
    save(OutputFile, 'onsets', 'offsets', 'labels', 'min_int', 'min_dur', 'sm_win', 'threshold');
end

% Calculated durations (durs) and gaps and the duration of the entire file
% (len)
durs = offsets - onsets;
gaps = onsets(2:end) - offsets(1:end-1);
len = length(Song) * 1000/Fs;