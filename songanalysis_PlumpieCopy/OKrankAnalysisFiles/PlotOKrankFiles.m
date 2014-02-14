function [data, SamplingRate] = PlotOKrankFiles(DirectoryName, FileName, ChanNo)

data = [];
SamplingRate = [];

cd(DirectoryName);

[datafid, message] = fopen(FileName, 'r');

[recfid, message] = fopen([FileName, '.rec'], 'r');

disp(FileName);

if ((recfid) > 0)
    while (~feof(recfid))
        tline = fgetl(recfid);
        if (strfind(tline, 'ai_freq'))
            ColonIndex = find(tline == ':');
            SamplingRate = str2double(tline((ColonIndex + 1):end))
        end

        if (strfind(tline, 'n_ai_chan'))
            ColonIndex = find(tline == ':');
            NoOfChannels = str2double(tline((ColonIndex + 1):end))
        end

        if (strfind(tline, 'n_samples'))
            ColonIndex = find(tline == ':');
            NoOfSamples = str2double(tline((ColonIndex + 1):end))
            break;
        end
    end

    fclose(recfid);

    %figure;
    %set(gcf, 'Color', 'w');
    %axes('position',[0.15 0.35 0.75 0.6]);
    %hold on;


    if ((ChanNo) > 100)
        for i = 1:NoOfChannels,
            fseek(datafid, (i - 1) * 2, 'bof');
            [data, num_read] = fread(datafid, inf, 'uint16', (NoOfChannels - 1) * 2);
            data = (data - 32768) * 10/32768;
            if (num_read ~= NoOfSamples)
                disp(['No of samples does not match that of recfile: ',FileName]);
            end
            time = 0:1/SamplingRate:(length(data)/SamplingRate);
            time(end) = [];
     %       plot(time, (data + (i-1)*5));
        end
    else
        i = ChanNo;
        fseek(datafid, (i - 1) * 2, 'bof');
        [data, num_read] = fread(datafid, inf, 'uint16', (NoOfChannels - 1) * 2);
        data = (data - 32768) * 10/32768;
        if (num_read ~= NoOfSamples)
            disp(['No of samples does not match that of recfile: ',FileName]);
        end
        time = 0:1/SamplingRate:(length(data)/SamplingRate);
        time(end) = [];
    %    plot(time, (data + (i-1)*5));
    end
    %axis('tight');

    %axes('position',[0.15 0.05 0.75 0.15]);
    %hold on;

    % fseek(datafid, (0) * 2, 'bof');
    % [data, num_read] = fread(datafid, inf, 'uint16', (NoOfChannels - 1) * 2);
    % data = (data - 32768) * 10/32768;
    %disp(num_read);

    % if (length(data) > 0;
    %     
    %     Fs = SamplingRate;
    % 
    %     % Now using an 8 pole butterworth bandpass filter as default.
    %     [b,a]=butter(8,[300*2/Fs, 10000*2/Fs]);
    % 
    %     FiltSong=filtfilt(b, a, data);
    % 
    %     if length(data) ~= length(FiltSong) 
    %       disp(['warning! bandpass: input and output file lengths do not match!']);
    %     end
    % 
    %     nfft=round(Fs*8/1000);
    %     nfft = 2^nextpow2(nfft);
    %     spect_win = hanning(nfft);
    %     noverlap = round(0.9*length(spect_win)); %number of overlapping points       
    %     %now calculate spectrogram
    %     [spect, freq, time] = spectrogram(FiltSong,spect_win,noverlap,nfft,Fs,'yaxis');
    %     idx_spect=scale_spect(spect);  %calculate index array for spectrogram
    %     f_min = freq(1);
    %     f_max = freq(length(freq));
    %     freq_spect = [f_min, f_max];
    %     t_min = time(1);
    %     t_max = time(end);
    % 
    %     %adjust time axis for spectrogram offset (1/2 window duration in ms)
    %     t_min = t_min + 0.5*nfft/Fs;  
    %     t_max = t_max + 0.5*nfft/Fs;  
    % 
    %     time_spect = [t_min, t_max];                
    % 
    %     %cm = disp_idx_spect(idx_spect, time_spect, freq_spect, -55, ...
    %     %    0, 1.5, 'gray', 'classic');
    %     %axis([0 time(end) 300 10000]);
    %     %set(gca,'Visible','off');
    % end
end

if ((datafid) > 0)
    fclose(datafid);
end
   


