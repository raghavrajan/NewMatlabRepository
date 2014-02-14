function [] = PlotFiringRateOverTime(DirectoryName, BirdName, Date, ChanNo, StartTime, EndTime)

cd(DirectoryName);

FigureNo = 0;
PlotNo = 1;
RowNo = 5;
ColNo = 3;
TotalPlots = RowNo * ColNo;
TimeInterval = 3;
PlotGap = 0.05;
PlotHeight = 0.9/RowNo - PlotGap;
PlotWidth = 0.9/ColNo - PlotGap;

TotalPlotNo = 1;
PresentRowNo = RowNo;
PresentColNo = 1;

for i = StartTime:TimeInterval:EndTime,
    Files = dir([BirdName,'_',Date,num2str(i),'*']);
    if (length(Files) > 0)
        [data, Fs] = PlotOKrankFiles(DirectoryName,Files(1).name,ChanNo);
        if (length(data) > 0)
            time = 0:1/Fs:length(data)/Fs;
            time(end) = [];
            
            if (mod(TotalPlotNo, TotalPlots) == 1)
                FigureNo = FigureNo + 1;
                PresentRowNo = RowNo;
                PresentColNo = 1;
                PlotNo = 1;
            end
            
            if (PlotNo > ColNo)
                PlotNo = 1;
                PresentRowNo = PresentRowNo - 1;
                PresentColNo = 1;
            end
            
            figure(FigureNo);
            set(gcf,'Color','w');
            AllPlots(TotalPlotNo) = axes('position',[(0.1 + ((PresentColNo - 1) * (PlotWidth + PlotGap))),(0.1 + ((PresentRowNo - 1) * (PlotHeight + PlotGap))), PlotWidth, PlotHeight]);
            plot(time(1:Fs*2), data(1:Fs*2));
            TitleString = [Files(1).name((end-5):(end-4)), ':', Files(1).name((end-3):(end-2)),':', Files(1).name((end-1):(end))];
            title(TitleString,'FontWeight','bold','FontSize',14);
            axis tight;
            PlotNo = PlotNo + 1;
            PresentColNo = PresentColNo + 1;
            TotalPlotNo = TotalPlotNo + 1;
        end
    end
end

PresentFigNo = 0;
for i = 1:length(AllPlots),
    if (mod(i,TotalPlots) == 1)
        PresentFigNo = PresentFigNo + 1;
    end
    figure(PresentFigNo);
    axes(AllPlots(i));
    temp(i,:) = axis;
end    

FinalPlotAxis = max(temp);
TempMin = min(temp);
FinalPlotAxis(3) = TempMin(3) * 1.05;
FinalPlotAxis(4) = FinalPlotAxis(4) * 1.05;

PresentFigNo = 0;
for i = 1:length(AllPlots),
    if (mod(i,TotalPlots) == 1)
        PresentFigNo = PresentFigNo + 1;
    end
    figure(PresentFigNo);
    axes(AllPlots(i));
    axis(FinalPlotAxis);
end    

disp('Finished');
            