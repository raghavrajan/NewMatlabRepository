function [PitchShiftedSound] = PitchShiftSound(RawData, Fs, FreqStretch)

if (FreqStretch ~= 1)
    % WarpedMotifRawData = pvoc(RawData, FreqStretch); 
    % Using another algorithm to shift pitch - seems to work pretty well
    % It expresses changes in pitch in terms of cents
    % To calculate cents I used the formula
    % cent change = 1200 * log(f2/f1) / log(2)
    parameter.fsAudio = Fs;
    PitchShiftedSound = pitchShiftViaTSM(RawData, 1200 * log(FreqStretch) / log(2), parameter); 
    % WarpedMotifRawData = resample(WarpedMotifRawData, length(RawData), length(WarpedMotifRawData));
else
    PitchShiftedSound = RawData;
end
