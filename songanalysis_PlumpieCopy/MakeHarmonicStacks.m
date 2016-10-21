function [Output] = MakeHarmonicStacks(NumStacks, OnsetTimes, Freqs, Fs, Durs, TotalDur, varargin)

if (nargin > 6)
    PlotToggle = varargin{1};
else
    PlotToggle = 'on';
end

if ((length(OnsetTimes) ~= NumStacks))
    msgbox('You have not specified the onset times for all of the harmonic stacks');
end

if (~isscalar(Freqs))
    if ((length(Freqs) ~= NumStacks))
        msgbox('You have not specified the fundamental frequencies (FF in Hz) for all of the harmonic stacks. Either specify for each one individually, or just put one value corresponding to the FF for all stacks');
    end
else
    Freqs = ones(size(OnsetTimes))*Freqs;
end

if (~isscalar(Durs))
    if ((length(Durs) ~= NumStacks))
        msgbox('You have not specified the durations (in seconds) for all of the harmonic stacks. Either specify for each one individually, or just put one value corresponding to the duration for all stacks');
    end
else
    Durs = ones(size(OnsetTimes))*Durs;
end
    
for i = 1:NumStacks,
    t = (1:round(Durs(i)*Fs))/Fs;
    HarmonicIndex = 1;
    Stack{i} = zeros(size(t));
    while (HarmonicIndex*Freqs(i) <= 8000)
        s = sin(2*pi*HarmonicIndex*Freqs(i)*t);
        HarmonicIndex = HarmonicIndex + 1;
        r = sin(linspace(0,pi/2,round(length(t)/10)));
        r = [r, ones(1, length(t) - length(r)*2), fliplr(r)];
        s = s.*r;
        Stack{i} = Stack{i} + s;
    end
end

Output = [];
for i = 1:length(OnsetTimes),
    if (i == 1)
        GapDur = OnsetTimes(i);
    else
        GapDur = OnsetTimes(i) - (OnsetTimes(i-1) + Durs(i-1));
    end
    Output = [Output [zeros(1, (round(GapDur*Fs))) Stack{i}]];
end
EndGapDur = TotalDur - (OnsetTimes(end) + Durs(end));
Output = [Output [zeros(1, (round(EndGapDur*Fs)))]];
Output = Output + rand(size(Output))/1000;
Output = Output/(max(Output)*1.1);

if (strfind(PlotToggle, 'on'))
    PlotSpectrogram_SongVar(Output, Fs);
end