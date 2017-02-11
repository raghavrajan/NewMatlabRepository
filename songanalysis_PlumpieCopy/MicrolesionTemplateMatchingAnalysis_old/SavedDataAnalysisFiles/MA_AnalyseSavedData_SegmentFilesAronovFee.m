function [onsets, offsets] = MA_AnalyseSavedData_SegmentFilesAronovFee(SongBout, Fs, min_int, min_dur, NoteFile, BoutOnset, BoutOffset)

sm_win = 2.5;

[LogAmplitude] = ASSLCalculateLogAmplitudeAronovFee(SongBout, Fs, (1:1:length(SongBout))/Fs, 8, 0.9);
[threshold] = ASSLCalculateFisherThreshold(LogAmplitude);
[onsets, offsets] = ASSLSegmentDataAronovFee(LogAmplitude, Fs, min_int, min_dur, threshold);

if (~isempty(onsets))
    labels = repmat('0', 1, length(onsets));
else
    labels = [];
end

save(NoteFile, 'onsets', 'offsets', 'labels', 'threshold', 'min_dur', 'min_int', 'BoutOnset', 'BoutOffset', 'sm_win');
            


