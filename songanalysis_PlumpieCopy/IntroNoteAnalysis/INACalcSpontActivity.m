function [SpontActivity] = INACalcSpontActivity(MotifIntroNoteIndices, ASpikeData)

BoutNo = 1;
if (~isfield(MotifIntroNoteIndices, 'Onsets'))
    SpontActivity.Spikes = [];
    SpontActivity.VocalInterval = [];
    SpontActivity.NoofIntroNotes = [];
    return;
end
for i = 1:length(MotifIntroNoteIndices.NoofINs),
    if ((MotifIntroNoteIndices.NoofINs{i} == 0))
        continue;
    end
    if (isempty(MotifIntroNoteIndices.NoofINs{i}))
        continue;
    end
    TempSpikeIndices = find((ASpikeData >= (MotifIntroNoteIndices.Onsets{i}(1) - MotifIntroNoteIndices.VocalInterval{i})) & (ASpikeData <= MotifIntroNoteIndices.Onsets{i}(1)));
    SpontActivity.Spikes{BoutNo} = ASpikeData(TempSpikeIndices) - MotifIntroNoteIndices.Onsets{i}(1);
    SpontActivity.VocalInterval{BoutNo} = MotifIntroNoteIndices.VocalInterval{i};
    SpontActivity.NoofIntroNotes{BoutNo} = MotifIntroNoteIndices.NoofINs{i};
    BoutNo = BoutNo + 1;
end
disp('Calculated spontaneous activity');