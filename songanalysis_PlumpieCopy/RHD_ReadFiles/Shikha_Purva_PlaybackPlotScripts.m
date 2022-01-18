function [] = Shikha_Purva_PlaybackPlotScripts(SongFile)

[Time, SpikeData, SongData, Fs] = read_Intan_RHD2000_file(SongFile);

% Plot Spectrogram and spike data in panels of 10s length

FileLen = length(SongData)/Fs.amplifier_sample_rate;
[LogAmplitude] = ASSLCalculateLogAmplitudeAronovFee(SongData - mean(SongData), Fs.amplifier_sample_rate, Time, 8, 0.5);

figure;
p = panel();
p.pack(6, 2);

PlotDur = 11; % in seconds
StartIndex = 1;
PlotIndex = 0;
while (StartIndex < length(SongData))
    PlotIndex = PlotIndex + 1;
    if (PlotIndex < 7)
        p(PlotIndex, ceil(PlotIndex/6)).select();
    else
        p(PlotIndex-6, ceil(PlotIndex/6)).select();
    end
    EndIndex = StartIndex + (PlotDur*Fs.amplifier_sample_rate);
    if (EndIndex > length(SongData))
        EndIndex = length(SongData);
    end
    PlotSpectrogramInAxis_SongVar(SongData(StartIndex:EndIndex), Time(StartIndex:EndIndex), Fs.amplifier_sample_rate, gca);
    
    
    PlotIndex = PlotIndex + 1;
    if (PlotIndex < 7)
        p(PlotIndex, ceil(PlotIndex/6)).select();
    else
        p(PlotIndex-6, ceil(PlotIndex/6)).select();
    end
    [Ax, H1, H2] = plotyy(Time(StartIndex:EndIndex), SpikeData(StartIndex:EndIndex), Time(StartIndex:EndIndex), LogAmplitude(StartIndex:EndIndex));
    axes(Ax(1));
    axis tight;
    axes(Ax(2));
    axis tight;
    
    StartIndex = EndIndex;
end
    
