function [] = PlotSongSAPFeatures(DirectoryName,FileName,TimeLength,FeatureName)

ChannelNo = 0;

% Get Song Data from observer file
channel_string = strcat('obs',num2str(0),'r');
[rawsong,Fs] = soundin_copy(DirectoryName,FileName,channel_string);

% Convert to uV - 5V on the data acquisition is 32768
rawsong = rawsong * 1/32768;

% Time axis
time = 0:1/32000:(length(rawsong)/32000);
time(end) = [];

PlotWidth = 0.85;
PlotGap = 0.05;
PlotHeight = (0.85 - 3*PlotGap)/3;
index = 1;
PlotIndex = 1;

for i = 1:ceil(time(end)/TimeLength),
    
    if (mod(i-1,3) == 0)
        PlotFigure(index) = figure;
        set(gcf,'Color','w');
        index = index + 1;
        PlotIndex = 1;
    end

    
    indices = find((time > ((i-1)*TimeLength)) & (time < (i * TimeLength)));

    % Plot song spectrogram
    SongPlot(PlotIndex) = axes('position',[0.1 (0.1 + ((PlotIndex - 1) * (PlotGap + PlotHeight))) PlotWidth PlotHeight]);
    PlotIndex = PlotIndex + 1;
    set(SongPlot,'FontSize',16);

    filtsong = bandpass(rawsong(indices),Fs,300,10000);
    
    nfft=round(Fs*8/1000);
    nfft = 2^nextpow2(nfft);
    spect_win = hanning(nfft);
    noverlap = round(0.9*length(spect_win)); %number of overlapping points       
    %now calculate spectrogram
%     [spect, freq, time_song] = specgram(filtsong, nfft, Fs, spect_win, noverlap);
    [spect, freq, time_song] = spectrogram(filtsong,spect_win,noverlap,nfft,Fs,'yaxis');
    idx_spect=scale_spect(spect);  %calculate index array for spectrogram
    f_min = freq(1);
    f_max = freq(length(freq));
    freq_spect = [f_min, f_max];

    t_min = time(indices(1)); %convert to ms
    t_max = time(indices(end)); %convert to ms
    %adjust time axis for spectrogram offset (1/2 window duration in ms)
    %t_min = t_min + 0.5*nfft*1000/Fs;  
    %t_max = t_max + 0.5*nfft*1000/Fs;  
    t_min = t_min;  
    t_max = t_max;  
    time_spect = [t_min, t_max];                

    cm = disp_idx_spect(idx_spect, time_spect, freq_spect, -70, ...
        5, 1.5, 'gray', 'classic');
    axis([((i-1)*TimeLength) (i*TimeLength) 300 10000]);
    
    % ylabel('Frequency (Hz)','FontSize',16);
    
    [m_spec_deriv , m_AM, m_FM ,m_Entropy , m_amplitude , m_Freq, m_PitchGoodness , m_Pitch , Pitch_chose , Pitch_weight ]=segment_data(rawsong(indices),Fs);

    load('Parameters');

    Axes2(PlotIndex - 1) = axes('position',get(SongPlot(PlotIndex - 1),'position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','r','YColor','r');
    
    if (strfind(FeatureName,'fm'))
        tmp_plot_function = m_FM/90*250*44100/1024;
        FeatureLine(PlotIndex -1) = line((1:size(m_spec_deriv,1))/1000,tmp_plot_function,'Color','r','Parent',Axes2(PlotIndex - 1));            
    else
        if (strfind(FeatureName,'pitch'))
            FeatureLine = line((1:size(m_spec_deriv,1))/1000,Pitch_chose,'Color','r','Parent',Axes2(PlotIndex - 1));            
        else
            if (strfind(FeatureName,'pgoodness'))
                tmp_plot_function = m_PitchGoodness - min(m_PitchGoodness);
                tmp_plot_function = tmp_plot_function/max(tmp_plot_function)*param.y_length* 32000/1024;
                FeatureLine = line((1:size(m_spec_deriv,1))/1000,tmp_plot_function,'Color','r','Parent',Axes2(PlotIndex - 1));            
            else
                if (strfind(FeatureName,'entropy'))
                    tmp_plot_function=m_Entropy - min(m_Entropy);
                    tmp_plot_function = tmp_plot_function/max(tmp_plot_function)*param.y_length*32000/1024;                
                    FeatureLine = line((1:size(m_spec_deriv,1))/1000,tmp_plot_function,'Color','r','Parent',Axes2(PlotIndex - 1));            
                end
            end
        end
    end
    axis([((i-1)*TimeLength) (i*TimeLength) min(tmp_plot_function) max(tmp_plot_function)]);
end

disp('Finished');

