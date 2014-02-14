function [] = PlotOKrankData_ReferenceSubtracted(DirectoryName, FileName, RefChan, varargin)

if (nargin > 3)
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
    SubPlotAxes = axes('position',[0.15 0.15 0.75 0.8]);
    hold on;

    DataMax = 0;

    fseek(datafid, (RefChan - 1) * 2, 'bof');
    [RefData, num_read] = fread(datafid, inf, 'uint16', (NoOfChannels - 1) * 2);
    RefData = (RefData - 32768) * 10/32768;
    
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
            plot_motif_spectrogram(data,Fs, MainFigure, SubPlotAxes)
            DataMax = 10000;
            data = data * 2000;
            data = data + DataMax + 1000 + abs(min(data));
            DataMax = max(data);
            figure(MainFigure);
            axes(SubPlotAxes);
            hold on;
            plot(time, data);
        else
            if (exist('ChanNo','var'))
                for (j = 1:length(ChanNo)),
                    if (ChanNo(j) == i)
                        data = data * 2000;
                        data = data + DataMax + 1000 + abs(min(data));
                        DataMax = max(data);
                        figure(MainFigure);
                        axes(SubPlotAxes);
                        hold on;
                        plot(time, data);
                    end
                end
            else
                data = data - RefData;
                data = data * 2000;
                data = data + DataMax + 1000 + abs(min(data));
                DataMax = max(data);
                figure(MainFigure);
                axes(SubPlotAxes);
                hold on;
                plot(time, data);
                axis tight;
            end
        end
    end
    axis tight;
    
end

if ((datafid) > 0)
    fclose(datafid);
end
   


