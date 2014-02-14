function [] = PlotUWEventParameters(DirFileInfo, UnDirFileInfo, MedianMotif, MainFigure)

figure(MainFigure);

axes('Position',[0.59 0.82 0.15 0.1]);
set(gca,'Box','off');
hold on;

MaxNoofSpikes = 0;

if ((isfield(DirFileInfo.UWEventParameters, 'NoofSpikes')) || (isfield(UnDirFileInfo.UWEventParameters, 'NoofSpikes')))
    if (isfield(DirFileInfo.UWEventParameters,'NoofSpikes'))
        if (length(DirFileInfo.UWEventParameters.NoofSpikes) > 0)
            for i = 1:size(DirFileInfo.UWEventParameters.NoofSpikes,2),
                errorbar(DirFileInfo.UWEventParameters.EventTimes(i),mean(DirFileInfo.UWEventParameters.NoofSpikes(:,i)),std(DirFileInfo.UWEventParameters.NoofSpikes(:,i)),'ks');
                MaxNoofSpikes = max(MaxNoofSpikes, (mean(DirFileInfo.UWEventParameters.NoofSpikes(:,i)) + std(DirFileInfo.UWEventParameters.NoofSpikes(:,i))));
            end
        end
    end

    if (isfield(UnDirFileInfo.UWEventParameters,'NoofSpikes'))
        if (length(UnDirFileInfo.UWEventParameters.NoofSpikes) > 0)
            for i = 1:size(UnDirFileInfo.UWEventParameters.NoofSpikes,2),
                errorbar(UnDirFileInfo.UWEventParameters.EventTimes(i),mean(UnDirFileInfo.UWEventParameters.NoofSpikes(:,i)),std(UnDirFileInfo.UWEventParameters.NoofSpikes(:,i)),'Color',[0.6 0.6 0.6],'Marker','d','MarkerEdgeColor',[0.6 0.6 0.6]);
                NoofEvents = size(UnDirFileInfo.UWEventParameters.NoofSpikes,2);                
                MaxNoofSpikes = max(MaxNoofSpikes, (mean(UnDirFileInfo.UWEventParameters.NoofSpikes(:,i)) + std(UnDirFileInfo.UWEventParameters.NoofSpikes(:,i))));            
            end
        end
    end
    axis([0 (MedianMotif.Length) 0 (MaxNoofSpikes * 1.25)]);
end

set(gca,'FontSize',8,'FontWeight','bold');
xlabel('Event time (sec)','FontSize',8,'FontWeight','bold');
ylabel('No of Spikes','FontSize',8,'FontWeight','bold');    
title('No of spikes per event','FontSize',8,'FontWeight','bold');    

figure(MainFigure);

axes('Position',[0.81 0.82 0.15 0.1]);
set(gca,'Box','off');
hold on;

MaxJitter = 0;

if ((isfield(DirFileInfo.UWEventParameters, 'NoofSpikes')) || (isfield(UnDirFileInfo.UWEventParameters, 'NoofSpikes')))
    if (isfield(DirFileInfo.UWEventParameters,'NoofSpikes'))
        if (length(DirFileInfo.UWEventParameters.NoofSpikes) > 0)
            for i = 1:size(DirFileInfo.UWEventParameters.NoofSpikes,2),
                plot(DirFileInfo.UWEventParameters.EventTimes(i),DirFileInfo.UWEventParameters.JittertoMeanFirstSpikeTime(i)*1000,'ks');
                MaxJitter = max(MaxJitter,DirFileInfo.UWEventParameters.JittertoMeanFirstSpikeTime(i));
            end
        end
    end
    if (isfield(UnDirFileInfo.UWEventParameters,'NoofSpikes'))
        if (length(UnDirFileInfo.UWEventParameters.NoofSpikes) > 0)
            for i = 1:size(UnDirFileInfo.UWEventParameters.NoofSpikes,2),
                plot(UnDirFileInfo.UWEventParameters.EventTimes(i), UnDirFileInfo.UWEventParameters.JittertoMeanFirstSpikeTime(i)*1000,'Marker','d','MarkerEdgeColor',[0.6 0.6 0.6]);
                MaxJitter = max(MaxJitter,UnDirFileInfo.UWEventParameters.JittertoMeanFirstSpikeTime(i));                
            end
        end
    end
    axis([0 (MedianMotif.Length) 0 (MaxJitter * 1000 * 1.25)]);
end

set(gca,'FontSize',8,'FontWeight','bold');
xlabel('Event Time (sec)','FontSize',8,'FontWeight','bold');
ylabel('Jitter (ms)','FontSize',8,'FontWeight','bold');    
title(['  Jitter of first spike '; 'to mean first spike time'],'FontSize',8,'FontWeight','bold');