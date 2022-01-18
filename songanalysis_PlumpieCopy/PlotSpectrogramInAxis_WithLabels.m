function [] = PlotSpectrogramInAxis_WithLabels(pathname, filename, FileType, gca, ColourMap, varargin)

PresentDir = pwd;

Slash = find((filename == '/') | (filename == '\'));
if (~isempty(Slash))
    filename = filename(Slash(end)+1:end);
end

if (nargin > 5),
    TimeRange = varargin{1};
end

if (nargin > 6),
    LabelFontSize = varargin{2};
else
    LabelFontSize = 8;
end

if (nargin > 7)
    LabelFontWeight = varargin{3};
else
    LabelFontWeight = 'normal';
end

if (nargin > 8)
    PlotScaleBar = varargin{4};
else
    PlotScaleBar = 0;
end

ChannelNo = 1;

% Get Song Data from observer file

if (strfind(FileType,'obs'))
    channel_string = strcat('obs',num2str(0),'r');
    ChannelNo = 0;
else
    channel_string = ChannelNo;
end

[rawsong, Fs] = GetData(pathname, filename, FileType, ChannelNo);

% Time axis
time = 0:1/Fs:(length(rawsong)/Fs);
time(end) = [];

if (exist('TimeRange', 'var'))
    start_time_index = find(time <= TimeRange(1), 1, 'last');
    end_time_index = find(time <= TimeRange(2), 1, 'last');
else
    start_time_index = 1;
    end_time_index = length(time);
end

% Plot song spectrogram

time = time(start_time_index:end_time_index);
rawsong = rawsong(start_time_index:end_time_index);

filtsong = bandpass(rawsong,Fs,300,10000);
%filtsong = bandpass_fft_filter(rawsong,8000,1500,Fs);
%filtnoise = bandpass_fft_filter(rawsong,1200,100,Fs);
Len = round(Fs*2/1000);
h = ones(1, Len)/Len;

%SquaredSong = FindEnvelope(time, filtsong, filtsong);
%SmoothSong = conv(h, SquaredSong);
%Offset = round((length(SmoothSong) - length(filtsong))/2);
%SmoothSong = SmoothSong((1 + Offset):(length(filtsong) + Offset));

% SquaredNoise = FindEnvelope(time, filtnoise, filtnoise);
% SmoothNoise = conv(h, SquaredNoise);
% Offset = round((length(SmoothNoise) - length(filtnoise))/2);
% SmoothNoise = SmoothNoise((1 + Offset):(length(filtnoise) + Offset));
% 
% figure;
% plot(time,SmoothSong,'r');
% hold on;
% plot(time,SmoothNoise,'b');
% 
nfft=round(Fs*8/1000);
nfft = 2^nextpow2(nfft);
spect_win = hanning(nfft);
noverlap = round(0.95*length(spect_win)); %number of overlapping points       

%now calculate spectrogram
%     [spect, freq, time_song] = specgram(filtsong, nfft, Fs, spect_win, noverlap);
[spect, freq, time_song] = spectrogram(filtsong,spect_win,noverlap,nfft,Fs,'yaxis');
idx_spect=scale_spect(spect);  %calculate index array for spectrogram
f_min = freq(1);
f_max = freq(length(freq));
freq_spect = [f_min, f_max];

t_min = time(1); %convert to ms
t_max = time(end); %convert to ms
%adjust time axis for spectrogram offset (1/2 window duration in ms)
%t_min = t_min + 0.5*nfft/Fs;  
%t_max = t_max + 0.5*nfft/Fs;  
time_spect = [t_min, t_max];   
axes(gca);
hold off;
cm = disp_idx_spect(idx_spect, time_spect, freq_spect, -55, ...
        0, 2, ColourMap, 'classic');
axis([t_min t_max 300 8000]);
set(gca, 'FontSize', 10);
% set(gca, 'XTick', []);
xlabel('Time (sec)', 'FontSize', 12);
ylabel('Frequency (Hz)', 'FontSize', 12);
zoom xon;
hold on;
%plot(time, (filtsong * 1000) + 12000);
%plot(time, (filtsong2 * 1000) + 14000,'r');
%NoteTimes = SmoothSong > 10*RMS;
%plot(time,NoteTimes * 14000,'k');

Notes = load(fullfile(pathname, 'ASSLNoteFiles', [filename, '.not.mat']));
for i = 1:length(Notes.onsets),
    % plot([Notes.onsets(i)/1000 Notes.onsets(i)/1000 Notes.offsets(i)/1000 Notes.offsets(i)/1000], [300 8000 8000 300], 'b', 'LineWidth', 2);
    if (nargin > 5)
        if ((Notes.onsets(i)/1000 >= TimeRange(1)) && (Notes.offsets(i)/1000 <= TimeRange(2)))
            text((mean([Notes.onsets(i) Notes.offsets(i)]))/1000, 7500, Notes.labels(i), 'FontSize', LabelFontSize, 'Color', 'b', 'FontWeight', LabelFontWeight, 'HorizontalAlignment', 'center');
        end
    else
        text((mean([Notes.onsets(i) Notes.offsets(i)]))/1000, 7500, Notes.labels(i), 'FontSize', LabelFontSize, 'Color', 'b', 'FontWeight', LabelFontWeight, 'HorizontalAlignment', 'center');
    end
end

if (PlotScaleBar == 1)
    plot([(t_min+0.1) (t_min+0.2)], [300 300], 'k', 'LineWidth', 1);
    text((t_min + 0.1), -200, '100 ms', 'FontSize', 6, 'HorizontalAlignment', 'left');
end
set(gcf, 'Color', 'w');
set(gca, 'XColor', 'w');
set(gca, 'YColor', 'w');
%disp('Finished');