function [] = plot_observer_data(pathname,filename,ChannelNo,varargin);

if (length(varargin) > 0)
    start_time = varargin{1};
    end_time = varargin{2};
end

% Get Spike Data from observer file
channel_string = strcat('obs',num2str(ChannelNo),'r');
[rawdata,Fs] = soundin_copy(pathname,filename,channel_string);

% Convert to uV - 5V on the data acquisition is 32768 and 10000 is the gain
rawdata = rawdata * 500/32768;

% Time axis
time = 0:1/32000:(length(rawdata)/32000);
time(end) = [];

% Get Song Data from observer file
channel_string = strcat('obs',num2str(0),'r');
[rawsong,Fs] = soundin_copy(pathname,filename,channel_string);

% Convert to uV - 5V on the data acquisition is 32768
rawsong = rawsong * 1/32768;

PlotFigure = figure;
set(gcf,'Color','w');

% Plot spike data
SpikeDataPlot = axes('position',[0.15 0.1 0.8 0.3]);
set(SpikeDataPlot,'FontSize',16);
plot(time,rawdata,'k');
xlabel('Time (sec)','FontSize',16);
ylabel('Amplitude (\muV)','FontSize',16);
if (exist('start_time','var'))
    indices = (find((time > start_time) & (time < end_time)));
    axis([start_time end_time min(rawdata(indices)) max(rawdata(indices))]);
else
    axis([0 time(end) min(rawdata) max(rawdata)]);
end

% Plot song amplitude data
SongAmplitudePlot = axes('position',[0.15 0.425 0.8 0.25]);
set(SongAmplitudePlot,'FontSize',16);
plot(time,rawsong,'k');
% ylabel('Amplitude (V)','FontSize',16);
if (exist('start_time','var'))
    axis([start_time end_time min(rawsong(indices)) max(rawsong(indices))]);
else
    axis([0 time(end) min(rawsong) max(rawsong)]);
end
set(gca,'xtick',[]);
% set(gca,'ytick',[]);

% Plot song spectrogram
SongSpectrogramPlot = axes('position',[0.15 0.7 0.8 0.25]);
set(SongSpectrogramPlot,'FontSize',16);

if (exist('start_time','var'))
    filtsong = bandpass(rawsong(indices),Fs,300,10000);
else
    filtsong = bandpass(rawsong,Fs,300,10000);
end

nfft=round(Fs*8/1000);
nfft = 2^nextpow2(nfft);
spect_win = hanning(nfft);
noverlap = round(0.9*length(spect_win)); %number of overlapping points       
%now calculate spectrogram
[spect, freq, time_song] = specgram(filtsong, nfft, Fs, spect_win, noverlap);
idx_spect=scale_spect(spect);  %calculate index array for spectrogram
f_min = freq(1);
f_max = freq(length(freq));
freq_spect = [f_min, f_max];

if (exist('start_time','var'))
    time = time(indices);
else
    time = time;
end
t_min = time(1); %convert to ms
t_max = time(length(time)); %convert to ms
%adjust time axis for spectrogram offset (1/2 window duration in ms)
%t_min = t_min + 0.5*nfft*1000/Fs;  
%t_max = t_max + 0.5*nfft*1000/Fs;  
t_min = t_min;  
t_max = t_max;  
time_spect = [t_min, t_max];                

cm = disp_idx_spect(idx_spect, time_spect, freq_spect, -55, ...
    5, 1.5, 'gray', 'classic');
axis([t_min t_max 300 10000]);
set(gca,'ytick',[]);
set(gca,'Visible','off');
set(gca,'xtick',[]);
% ylabel('Frequency (Hz)','FontSize',16);

disp('Finished');