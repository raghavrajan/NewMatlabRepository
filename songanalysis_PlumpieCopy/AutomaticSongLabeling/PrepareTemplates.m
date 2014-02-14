function [] = PrepareTemplates(DirectoryName, SongFileList, FileType)

if (DirectoryName(end) ~= '/')
    DirectoryName(end + 1) = '/';
end

cd(DirectoryName);

fid = fopen(SongFileList);

NoteNo = 0;

while (1)
    SongFileName = fgetl(fid);
    if (~ischar(SongFileName))
        break;
    end
    
    if (strfind(FileType, 'okrank'))
        [RawData, Fs] = ReadOKrankData(DirectoryName, SongFileName, 0);
    else
        if (strfind(FileType, 'obs'))
            [RawData, Fs] = soundin_copy(DirectoryName, SongFileName, 'obs0r');    
            RawData = RawData/max(RawData);
            RawData = RawData - mean(RawData);
        end
    end
    
    Time = (1:1:length(RawData))/Fs;
    
    WindowIncrement = 128;
    WindowLength = 256;
    WindowIndices = zeros(WindowLength, (ceil((length(RawData) - WindowLength)/WindowIncrement) + 1));
    WindowIndices(1,1) = 1;
    for i = 2:size(WindowIndices,2),
        WindowIndices(1,i) = WindowIndices(1,i-1) + WindowIncrement;
    end
    for i = 2:size(WindowIndices,1),
        WindowIndices(i,:) = WindowIndices(i-1,:) + 1;
    end
    
    RawData(end+1:WindowIndices(end,end)) = 0;
    FFTData = fft(RawData(WindowIndices));
    
    Freq = Fs/2*linspace(0, 1, WindowLength/2+1);
    FreqRows = find((Freq >=1700) & (Freq <= 7100));
    Amplitude = sum(log(2*abs(FFTData(FreqRows, :))));
    New_Fs = 1/(Time(WindowIndices(1,2)) - Time(WindowIndices(1,1)));
    
    XGauss = -32:1:32;
    STDGaussWin = 0.01 * New_Fs;
    GaussWin = 1/((STDGaussWin) * sqrt(2*pi)) * exp(-(XGauss.*XGauss)/(2 * STDGaussWin * STDGaussWin));
    SmoothAmplitude = conv(GaussWin, Amplitude);
    SmoothAmplitude = SmoothAmplitude(33:(length(SmoothAmplitude) - 32));
    
    figure;
    AmplitudePlot = axes('Position', [0.1 0.1 0.8 0.35]);
    hold on;
    plot(Time, RawData);
    plot(Time(WindowIndices(1,:)), Amplitude/max(Amplitude), 'r');
    plot(Time(WindowIndices(1,:)), SmoothAmplitude/max(SmoothAmplitude), 'g');
    plot(Time(WindowIndices(1,1:end-1)), diff(Amplitude)/max(diff(Amplitude)), 'k');
    axis tight;
    SpectrogramPlot = axes('Position', [0.1 0.7 0.8 0.25]);
    plot_motif_spectrogram(RawData, Fs, gcf, SpectrogramPlot);
    LabelPlot = axes('Position', [0.1 0.55 0.8 0.05]);
end
    
disp('Finished');