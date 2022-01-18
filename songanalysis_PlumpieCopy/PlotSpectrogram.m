function [MaxSpectVal] = PlotSpectrogram(pathname,filename,FileType,ColourMap, varargin)


Slash = find((filename == '/') | (filename == '\'));
if (~isempty(Slash))
    filename = filename(Slash(end)+1:end);
end

ChannelNo = 0;

% Get Song Data from observer file
if (strfind(FileType,'obs'))
    channel_string = strcat('obs',num2str(0),'r');
    [rawsong,Fs] = soundin_copy(pathname,filename,channel_string);
    
    % Convert to uV - 5V on the data acquisition is 32768
    rawsong = rawsong * 5/32768;
else
    if (strfind(FileType,'wav'))
        PresentDir = pwd;
        cd(pathname);
        [rawsong, Fs] = audioread(filename);
        cd(PresentDir);
    else 
        if (strfind(FileType, 'okrank'))
            [rawsong, Fs] = ReadOKrankData(pathname, filename, 0);
        end
    end
end


% Time axis
time = 0:1/Fs:(length(rawsong)/Fs);
time(end) = [];

if (length(varargin) > 0)
    start_time = varargin{1};
    end_time = varargin{2};
    start_time_index = find(time <= start_time,1,'last');
    end_time_index = find(time <= end_time,1,'last');
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
%[PowSpect, freq, spect, time_song] = CalculateMultiTaperSpectrogram(filtsong, Fs, 8, 1, 1.5);
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
figure;
set(gcf, 'Color', 'w');
set(gcf, 'Position', [153 440 1200 250]);
axes('Position',[0.08 0.25 0.88 0.6]);
cm = disp_idx_spect(idx_spect, time_spect, freq_spect, -55, ...
        -10, 1.2, ColourMap, 'classic');
axis([t_min t_max 300 8000]);
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Time (sec)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Frequency (kHz)', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'YTick', [0:2000:8000], 'YTickLabel', [0:2:8]);
title(filename, 'FontSize', 14, 'FontWeight', 'bold');
zoom xon;
hold on;
%plot(time, (filtsong * 1000) + 12000);
%plot(time, (filtsong2 * 1000) + 14000,'r');
%NoteTimes = SmoothSong > 10*RMS;
%plot(time,NoteTimes * 14000,'k');
disp('Finished');