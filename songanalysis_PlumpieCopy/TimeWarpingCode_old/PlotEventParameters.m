function [] = PlotEventParameters(DirFileInfo, UnDirFileInfo, MedianMotif, MainFigure)

figure(MainFigure);

axes('Position',[0.59 0.82 0.15 0.1]);
set(gca,'Box','off');
hold on;

MaxNoofSpikes = 0;

if ((isfield(DirFileInfo.EventParameters, 'NoofSpikes')) || (isfield(UnDirFileInfo.EventParameters, 'NoofSpikes')))
    if (isfield(DirFileInfo.EventParameters,'NoofSpikes'))
        if (length(DirFileInfo.EventParameters.NoofSpikes) > 0)
            for i = 1:size(DirFileInfo.EventParameters.NoofSpikes,2),
                errorbar(DirFileInfo.EventParameters.EventTimes(i),mean(DirFileInfo.EventParameters.NoofSpikes(:,i)),std(DirFileInfo.EventParameters.NoofSpikes(:,i)),'ks');
                MaxNoofSpikes = max(MaxNoofSpikes, (mean(DirFileInfo.EventParameters.NoofSpikes(:,i)) + std(DirFileInfo.EventParameters.NoofSpikes(:,i))));
            end
        end
    end

    if (isfield(UnDirFileInfo.EventParameters,'NoofSpikes'))
        if (length(UnDirFileInfo.EventParameters.NoofSpikes) > 0)
            for i = 1:size(UnDirFileInfo.EventParameters.NoofSpikes,2),
                errorbar(UnDirFileInfo.EventParameters.EventTimes(i),mean(UnDirFileInfo.EventParameters.NoofSpikes(:,i)),std(UnDirFileInfo.EventParameters.NoofSpikes(:,i)),'Color',[0.6 0.6 0.6],'Marker','d','MarkerEdgeColor',[0.6 0.6 0.6]);
                NoofEvents = size(UnDirFileInfo.EventParameters.NoofSpikes,2);                
                MaxNoofSpikes = max(MaxNoofSpikes, (mean(UnDirFileInfo.EventParameters.NoofSpikes(:,i)) + std(UnDirFileInfo.EventParameters.NoofSpikes(:,i))));            
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

if ((isfield(DirFileInfo.EventParameters, 'NoofSpikes')) || (isfield(UnDirFileInfo.EventParameters, 'NoofSpikes')))
    if (isfield(DirFileInfo.EventParameters,'NoofSpikes'))
        if (length(DirFileInfo.EventParameters.NoofSpikes) > 0)
            for i = 1:size(DirFileInfo.EventParameters.NoofSpikes,2),
                plot(DirFileInfo.EventParameters.EventTimes(i),DirFileInfo.EventParameters.JittertoMeanFirstSpikeTime(i)*1000,'ks');
                MaxJitter = max(MaxJitter,DirFileInfo.EventParameters.JittertoMeanFirstSpikeTime(i));
            end
        end
    end
    if (isfield(UnDirFileInfo.EventParameters,'NoofSpikes'))
        if (length(UnDirFileInfo.EventParameters.NoofSpikes) > 0)
            for i = 1:size(UnDirFileInfo.EventParameters.NoofSpikes,2),
                plot(UnDirFileInfo.EventParameters.EventTimes(i), UnDirFileInfo.EventParameters.JittertoMeanFirstSpikeTime(i)*1000,'Marker','d','MarkerEdgeColor',[0.6 0.6 0.6]);
                MaxJitter = max(MaxJitter,UnDirFileInfo.EventParameters.JittertoMeanFirstSpikeTime(i));                
            end
        end
    end
    axis([0 (MedianMotif.Length) 0 (MaxJitter * 1000 * 1.25)]);
end

set(gca,'FontSize',8,'FontWeight','bold');
xlabel('Event Time (sec)','FontSize',8,'FontWeight','bold');
ylabel('Jitter (ms)','FontSize',8,'FontWeight','bold');    
title(['  Jitter of first spike '; 'to mean first spike time'],'FontSize',8,'FontWeight','bold');