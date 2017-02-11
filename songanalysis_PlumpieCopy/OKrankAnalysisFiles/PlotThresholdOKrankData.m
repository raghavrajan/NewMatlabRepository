function [] = PlotThresholdOKrankData(DirectoryName, FileName, varargin)

if (nargin > 2)
    ChanNo = varargin{1};
end

cd(DirectoryName);

[datafid, message] = fopen(FileName, 'r');

[recfid, message] = fopen([FileName, '.rec'], 'r');

disp(FileName);

if ((recfid) > 0)
    while (~feof(recfid))
        tline = fgetl(recfid);
        if (strfind(tline, 'ai_freq'))
            ColonIndex = find(tline == ':');
            Fs = str2double(tline((ColonIndex + 1):end))
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

    MainFigure = figure;
    set(gcf, 'Color', 'w');
    SpectrogramAxes = axes('Position',[0.1 0.1 0.85 0.25]);
    SpikeAxes = axes('Position', [0.1 0.4 0.85 0.5]);
    hold on;

    DataMax = 0;

    for i = 1:NoOfChannels,
        fseek(datafid, (i - 1) * 2, 'bof');
        [data, num_read] = fread(datafid, inf, 'uint16', (NoOfChannels - 1) * 2);
        data = (data - 32768) * 10/32768;
        if (num_read ~= NoOfSamples)
            disp(['No of samples does not match that of recfile: ',FileName]);
        end
        time = 0:1/Fs:(length(data)/Fs);
        time(end) = [];
        
        if (i == 1)
            plot_motif_spectrogram(data,Fs, MainFigure, SpectrogramAxes)
            set(gca, 'FontSize', 14, 'FontWeight', 'bold')
            hold on;
        else
            if (exist('ChanNo','var'))
                for (j = 1:length(ChanNo)),
                    if (ChanNo(j) == i)
                        Prompt = {'Enter the lower threshold (uV)','Enter the upper threshold (uV)','Enter the window size for detecting spikes (no of points)','Enter the window size for skipping detection after a spike (no of points)'};
                        DialogTitle = 'Input for thresholding parameters';
                        DefaultParameters = {'-100','100','32','32'};

                        ThresholdingParameters = inputdlg(Prompt, DialogTitle, 1, DefaultParameters);

                        LowerThreshold = str2double(ThresholdingParameters{1});
                        UpperThreshold = str2double(ThresholdingParameters{2});
                        WindowSize = str2double(ThresholdingParameters{3});
                        SkipWindowSize = str2double(ThresholdingParameters{4});
                        
                        data = data * 100;

                        [TempSpikes, TempSpikeWaveforms] = FindSpikes(data,LowerThreshold,UpperThreshold,WindowSize,SkipWindowSize,0,Fs);
                        
                        figure(MainFigure);
                        axes(SpikeAxes);
                        hold on;
                        plot(time, data);
                        
                        for i = 1:length(TempSpikes),
                            Index = find(time < TempSpikes(i), 1, 'last');
                            plot(time((Index - 8):(Index + 23)), data((Index - 8):(Index + 23)),'r');
                        end
                        axis tight;
                    end
                end
            end
        end
    end
    axis tight;

    figure(MainFigure);
    title(FileName);
    axes(SpikeAxes);
    set(gca, 'FontSize', 14, 'FontWeight', 'bold')
    axis tight;
end

if ((datafid) > 0)
    fclose(datafid);
end
   


