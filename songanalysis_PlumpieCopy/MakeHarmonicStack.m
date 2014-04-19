function [Stack] = MakeHarmonicStack(Freq, Fs, Dur, TotalDur)

t = (1:round(Dur*Fs))/Fs;
HarmonicIndex = 1;
Stack = zeros(size(t));
while (HarmonicIndex*Freq <= 8000)
    s = sin(2*pi*HarmonicIndex*Freq*t);
    HarmonicIndex = HarmonicIndex + 1;
    r = sin(linspace(0,pi/2,round(length(t)/10)));
    r = [r, ones(1, length(t) - length(r)*2), fliplr(r)];
    s = s.*r;
    Stack = Stack + s;
end
Stack = [zeros(1,(round(TotalDur*Fs) - length(t))/2) Stack zeros(1,(round(TotalDur*Fs) - length(t))/2)];
Stack = Stack + rand(size(Stack))/1000;
Stack = Stack/(max(Stack)*1.1);
wavwrite(Stack, Fs, 16, [num2str(Freq), 'Hz_stack.wav']);
PlotSpectrogram(pwd, [num2str(Freq), 'Hz_stack.wav'], 'wav', 'hot');