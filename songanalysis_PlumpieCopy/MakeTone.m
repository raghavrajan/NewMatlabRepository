function [] = MakeTone(Freq, Fs, Dur, TotalDur)

t = (1:round(Dur*Fs))/Fs;
s = sin(2*pi*Freq*t);
r = sin(linspace(0,pi/2,round(length(t)/10)));
r = [r, ones(1, length(t) - length(r)*2), fliplr(r)];
s = s.*r;
s = [zeros(1,(round(TotalDur*Fs) - length(t))/2) s zeros(1,(round(TotalDur*Fs) - length(t))/2)];
s = s + rand(size(s))/100;
s = s/2;
wavwrite(s, Fs, 16, [num2str(Freq), 'Hz_tone.wav']);
PlotSpectrogram(pwd, [num2str(Freq), 'Hz_tone.wav'], 'wav');