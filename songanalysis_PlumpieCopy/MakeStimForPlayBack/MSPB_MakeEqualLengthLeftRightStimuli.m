function [ActualStimulus] = MSPB_MakeEqualLengthLeftRightStimuli(Stimulus, Silence, SilenceFs, MaxStimLength)

Stimulus = Stimulus(:);
Silence = Silence(:);

FinalStimLen = round(MaxStimLength * SilenceFs);
LenDiff = FinalStimLen - length(Stimulus);

if (LenDiff > 0)
    if (LenDiff > length(Silence))
        ActualStimulus = [Stimulus; repmat(Silence, floor(LenDiff/length(Silence)), 1)];
        ExtraDiff = FinalStimLen - length(ActualStimulus);
        ActualStimulus = [ActualStimulus; Silence(1:ExtraDiff)];
    end
else
    ActualStimulus = Stimulus;
end

ActualStimulus = ActualStimulus/sqrt(mean(ActualStimulus.*ActualStimulus));
ActualStimulus = ActualStimulus/max(ActualStimulus);

