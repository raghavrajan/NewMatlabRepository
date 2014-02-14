function  [Spectrogram] = PlotMotifSpectrogram(DirectoryName,MedianMotif,MainFigure,SpectrogramPlot)

[RawSong, Fs] = soundin_copy(DirectoryName,[MedianMotif.FileName{1}],'obs0r');
Time = 0:1/Fs:length(RawSong)/Fs;
Time(end) = [];

StartIndex = find(Time < (MedianMotif.StartTime),1,'last');
EndIndex = find(Time < (MedianMotif.StartTime + MedianMotif.Length),1,'last');


% Now using an 8 pole butterworth bandpass filter as default.
[b,a]=butter(8,[300*2/Fs, 10000*2/Fs]);

FiltSong=filtfilt(b, a, RawSong(StartIndex:EndIndex));
  
if length(RawSong(StartIndex:EndIndex)) ~= length(FiltSong) 
  disp(['warning! bandpass: input and output file lengths do not match!']);
end


nfft=round(Fs*8/1000);
nfft = 2^nextpow2(nfft);
spect_win = hanning(nfft);
noverlap = round(0.9*length(spect_win)); %number of overlapping points       
%now calculate spectrogram
[spect, freq, time] = spectrogram(FiltSong,spect_win,noverlap,nfft,Fs,'yaxis');
idx_spect=scale_spect(spect);  %calculate index array for spectrogram
f_min = freq(1);
f_max = freq(length(freq));
freq_spect = [f_min, f_max];
t_min = time(1);
t_max = time(end);

%adjust time axis for spectrogram offset (1/2 window duration in ms)
t_min = t_min + 0.5*nfft/Fs;  
t_max = t_max + 0.5*nfft/Fs;  

time_spect = [t_min, t_max];                

figure(MainFigure);
axes(SpectrogramPlot);
cm = disp_idx_spect(idx_spect, time_spect, freq_spect, -55, ...
    0, 2, 'gray', 'classic');
axis([-0.05 (MedianMotif.Length) 300 10000]);
set(gca,'Visible','off');
    