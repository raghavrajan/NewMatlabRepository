function [] = PlotFiringRateAtEachTimePoint(DirectoryName, BirdName, Date, ChanNo, StartTime, EndTime, LowerThreshold, UpperThreshold, WindowSize, SkipWindowSize)

cd(DirectoryName);

Files = dir([BirdName,'_',Date,'*']);
ActualFileNos = [];

for i = 1:length(Files),
    if (~length(strfind(Files(i).name, '.rec')))
        ActualFileNos(end + 1) = i;
    end
end

Files = Files(ActualFileNos);

for i = 1:length(Files),
    FileTimes(i) = str2double(Files(i).name((end - 5):end));
end

FileIndicesToAnalyse = find((FileTimes > StartTime) & (FileTimes < EndTime));
FilesToAnalyse = Files(FileIndicesToAnalyse);

FiringRate = [];
for i = 1:length(FilesToAnalyse),
    [RawData, Fs] = ReadOKrankData(DirectoryName, FilesToAnalyse(i).name, ChanNo);
    if (length(RawData) > 0)
        time = 0:1/Fs:length(RawData)/Fs;
        time(end) = [];
        TempSpikes = [];
        TempSpikeWaveforms = [];
        [TempSpikes, TempSpikeWaveforms] = FindSpikes(RawData(1:end),LowerThreshold,UpperThreshold,WindowSize,SkipWindowSize,0,Fs);
        FiringRate(end+1,1) = str2double(FilesToAnalyse(i).name((end-5):(end-4))) + str2double(FilesToAnalyse(i).name((end-3):(end-2)))/60 + str2double(FilesToAnalyse(i).name((end-1):(end)))/3600;
        FiringRate(end,2) = length(TempSpikes)/(time(end));
    end
end
    
figure;
set(gcf,'Color','w');
plot(FiringRate(:,1),FiringRate(:,2));
axis tight;
set(gca,'FontWeight','bold','FontSize',14);
xlabel('Time of day (hours)','FontWeight','bold','FontSize',14);
ylabel('Firing Rate (Hz)', 'FontWeight','bold','FontSize',14);
title([BirdName, ' ', Date],'FontWeight','bold','FontSize',14);

temp = axis;
hold on;
DrugTimes = inputdlg([{'Drug ON time'} {'Drug OFF time'}],'Drug ON, OFF times');

ColonIndex = find(DrugTimes{1} == ':');
DrugONTime = str2double(DrugTimes{1}(1:(ColonIndex - 1)));
Minutes = str2double(DrugTimes{1}((ColonIndex + 1):end))/60;
DrugONTime = DrugONTime + Minutes;
plot([DrugONTime DrugONTime],[0 temp(4)],'r-.');

ColonIndex = find(DrugTimes{2} == ':');
DrugOFFTime = str2double(DrugTimes{2}(1:(ColonIndex - 1)));
Minutes = str2double(DrugTimes{2}((ColonIndex + 1):end))/60;
DrugOFFTime = DrugOFFTime + Minutes;
plot([DrugOFFTime DrugOFFTime],[0 temp(4)],'r-.');

axis tight;
temp = axis;
temp(3) = 0;
axis(temp);

annotation('textbox',[0.05 0.05 0.8 0.0025],'FontSize',8,'FontWeight','bold','String',['Lower Threshold is ', num2str(LowerThreshold),' and Upper Threshold is ', num2str(UpperThreshold)], 'LineStyle','none');

disp('Finished');