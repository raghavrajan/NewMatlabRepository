function [] = PlotDIRUNDIR(BirdName,Directories,Motif);

for DirectoryNo = 1:size(Directories,1),
    cd(Directories(DirectoryNo,:));
    MotifDetails = [];
    for i = 1:length(Motif),
        FileName = [BirdName,'_',Directories(DirectoryNo,:),'_all_',Motif(i),'.mat'];
        MotifDetails(i).syllable = load(FileName);
    end
    cd ../;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The different columns of the cell array ffreq are as follows
    %   Column No.      Value
    %       1       Sound file
    %       2       Fundamental Frequency (using auto-correlation)
    %       3       mean Pitch (SAP measure - mean of the Pitch chose)
    %       4       mean Pitch Goodness (SAP measure)
    %       5       mean Frequency modulation Goodness (SAP measure)
    %       6       mean Entropy (SAP measure)
    %       7       mean Amplitude (SAP measure)
    %       8       mean Amplitude modulation (SAP measure)
    %       9       Syllable duration
    %      10       Time of day of the syllable
    %      11       Directed (D) or Undirected (U)
    %      12       The actual note file
    %      13       The bout no. of the motif in the file
    %      14       The motif no. within the bout that this motif is part of

    PlotTitles{1} = 'FF';
    PlotTitles{2} = 'Pitch';
    PlotTitles{3} = 'Pitch Goodness';
    PlotTitles{4} = 'FM';
    PlotTitles{5} = 'Entropy';
    PlotTitles{6} = 'Amplitude';
    PlotTitles{7} = 'AM';
    PlotTitles{8} = 'Duration';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now to do the plots
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PlotGap = 0.06;
    PlotWidth = 0.4;
    PlotHeight = (0.75 - 3*(PlotGap))/4;

    for Syllable = 1:length(Motif),
        if (exist('MainFeatureFigure','var'))
            if (length(MainFeatureFigure) >= Syllable)
            else
                MainFeatureFigure(Syllable) = figure;
                set(gcf,'Color','w');
                MainFigureTitle(Syllable) = annotation('textbox',[0.3 0.925 0.4 0.05]);
                set(MainFigureTitle(Syllable),'String',['Syllable ',Motif(Syllable)],'HorizontalAlignment','center','VerticalAlignment','middle','LineStyle','none','FontSize',16);
            end
        else
                MainFeatureFigure(Syllable) = figure;
                set(gcf,'Color','w');
                MainFigureTitle(Syllable) = annotation('textbox',[0.3 0.925 0.4 0.05]);
                set(MainFigureTitle(Syllable),'String',['Syllable ',Motif(Syllable)],'HorizontalAlignment','center','VerticalAlignment','middle','LineStyle','none','FontSize',16);
        end

        if (exist('FeatureTimeFigure','var'))
            if (length(FeatureTimeFigure) >= Syllable)
            else
                FeatureTimeFigure(Syllable) = figure;
                set(gcf,'Color','w');
                FeatureTimeFigureTitle(Syllable) = annotation('textbox',[0.3 0.925 0.4 0.05]);
                set(FeatureTimeFigureTitle(Syllable),'String',['Syllable ',Motif(Syllable)],'HorizontalAlignment','center','VerticalAlignment','middle','LineStyle','none','FontSize',16);
            end
        else
            FeatureTimeFigure(Syllable) = figure;
            set(gcf,'Color','w');
            FeatureTimeFigureTitle(Syllable) = annotation('textbox',[0.3 0.925 0.4 0.05]);
            set(FeatureTimeFigureTitle(Syllable),'String',['Syllable ',Motif(Syllable)],'HorizontalAlignment','center','VerticalAlignment','middle','LineStyle','none','FontSize',16);
        end

        if (exist('BoutNoFigure','var'))
            if (length(BoutNoFigure) >= Syllable)
            else    
                BoutNoFigure(Syllable) = figure;
                set(gcf,'Color','w');
                BoutNoFigureTitle(Syllable) = annotation('textbox',[0.3 0.925 0.4 0.05]);
                set(BoutNoFigureTitle(Syllable),'String',['Syllable ',Motif(Syllable)],'HorizontalAlignment','center','VerticalAlignment','middle','LineStyle','none','FontSize',16);
            end
        else
            BoutNoFigure(Syllable) = figure;
            set(gcf,'Color','w');
            BoutNoFigureTitle(Syllable) = annotation('textbox',[0.3 0.925 0.4 0.05]);
            set(BoutNoFigureTitle(Syllable),'String',['Syllable ',Motif(Syllable)],'HorizontalAlignment','center','VerticalAlignment','middle','LineStyle','none','FontSize',16);
        end

        if (exist('MotifNoFigure','var'))
            if (length(MotifNoFigure) >= Syllable)
            else    
                MotifNoFigure(Syllable) = figure;
                set(gcf,'Color','w');
                MotifNoFigureTitle(Syllable) = annotation('textbox',[0.3 0.925 0.4 0.05]);
                set(MotifNoFigureTitle(Syllable),'String',['Syllable ',Motif(Syllable)],'HorizontalAlignment','center','VerticalAlignment','middle','LineStyle','none','FontSize',16);
            end
        else
            MotifNoFigure(Syllable) = figure;
            set(gcf,'Color','w');
            MotifNoFigureTitle(Syllable) = annotation('textbox',[0.3 0.925 0.4 0.05]);
            set(MotifNoFigureTitle(Syllable),'String',['Syllable ',Motif(Syllable)],'HorizontalAlignment','center','VerticalAlignment','middle','LineStyle','none','FontSize',16);
        end

        DirectedIndices = find([MotifDetails(Syllable).syllable.ffreq{:,11}] == 'D');
        UnDirectedIndices = find([MotifDetails(Syllable).syllable.ffreq{:,11}] == 'U');    

        for i = 1:1,
           if (DirectoryNo == 1)
               if (mod(i,2) == 0)
                   figure(MainFeatureFigure(Syllable));
                   FeatureFigures(Syllable).Figure(i) = axes('position',[(0.1) (0.1 + ((i/2 - 1)*(PlotHeight + PlotGap))) PlotWidth PlotHeight]);
                   figure(FeatureTimeFigure(Syllable));
                   FeatureTimeFigures(Syllable).Figure(i) = axes('position',[(0.1) (0.1 + ((i/2 - 1)*(PlotHeight + PlotGap))) PlotWidth PlotHeight]);
                   figure(BoutNoFigure(Syllable));
                   BoutNoFigures(Syllable).Figure(i) = axes('position',[(0.1) (0.1 + ((i/2 - 1)*(PlotHeight + PlotGap))) PlotWidth PlotHeight]);
                   figure(MotifNoFigure(Syllable));
                   MotifNoFigures(Syllable).Figure(i) = axes('position',[(0.1) (0.1 + ((i/2 - 1)*(PlotHeight + PlotGap))) PlotWidth PlotHeight]);
               else
                   figure(MainFeatureFigure(Syllable));
                   FeatureFigures(Syllable).Figure(i) = axes('position',[(0.1 + PlotWidth + PlotGap) (0.1 + (((i - 1)/2)*(PlotHeight + PlotGap))) PlotWidth PlotHeight]);
                   figure(FeatureTimeFigure(Syllable));
                   FeatureTimeFigures(Syllable).Figure(i) = axes('position',[(0.1 + PlotWidth + PlotGap) (0.1 + (((i - 1)/2)*(PlotHeight + PlotGap))) PlotWidth PlotHeight]);
                   figure(BoutNoFigure(Syllable));
                   BoutNoFigures(Syllable).Figure(i) = axes('position',[(0.1 + PlotWidth + PlotGap) (0.1 + (((i - 1)/2)*(PlotHeight + PlotGap))) PlotWidth PlotHeight]);
                   figure(MotifNoFigure(Syllable));
                   MotifNoFigures(Syllable).Figure(i) = axes('position',[(0.1 + PlotWidth + PlotGap) (0.1 + (((i - 1)/2)*(PlotHeight + PlotGap))) PlotWidth PlotHeight]);
               end
           end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the features
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           axes(FeatureFigures(Syllable).Figure(i));
           set(gca,'FontSize',14);
           hold on;
           DirXValues = ones(length(DirectedIndices),1) + 2*(DirectoryNo - 1);
           UnDirXValues = 2 * ones(length(UnDirectedIndices),1) + 2*(DirectoryNo - 1);
           plot(DirXValues,[MotifDetails(Syllable).syllable.ffreq{DirectedIndices,(i+1)}],'r.');
           plot(UnDirXValues,[MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,(i+1)}],'b.');
           YMax = max(max([MotifDetails(Syllable).syllable.ffreq{DirectedIndices,(i+1)}]),max([MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,(i+1)}]));
           YMin = min(min([MotifDetails(Syllable).syllable.ffreq{DirectedIndices,(i+1)}]),min([MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,(i+1)}]));

           YRange = YMax - YMin;
           YMax = YMax + 0.1*YRange;
           YMin = YMin - 0.1*YRange;
           YRange = YMax - YMin;

           XMin = 0.5;
           XMax = 2 + 2*(DirectoryNo - 1) + 0.5;
           
           axis([XMin XMax YMin YMax]);

           if (i > 2)
               set(gca,'xtick',[]);
           else
               Xticks = 1:1:(2*DirectoryNo);
               for ticks = 1:length(Xticks),
                   if (mod(ticks,2) == 0)
                       Xticklabels{ticks} = 'Undirected';
                   else
                       Xticklabels{ticks} = 'Directed';
                   end
               end
               
               set(gca,'xtick',Xticks,'XTickLabel',Xticklabels,'FontSize',14);
               text((1 + 2*(DirectoryNo - 1)),(YMax),['(',num2str(length(DirectedIndices)),')']);
               text((2 + 2*(DirectoryNo - 1)),(YMax),['(',num2str(length(UnDirectedIndices)),')']);
           end
           axis auto;
           ReScaledAxis = ReScaleAxis(axis);
           axis(ReScaledAxis);
           TitleXCoord = ReScaledAxis(1) + (ReScaledAxis(2) - ReScaledAxis(1))/2;
           TitleYCoord = ReScaledAxis(4);
           title(PlotTitles{i},'FontSize',14,'Position',[TitleXCoord TitleYCoord]);
           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the features as a function of time of day

           axes(FeatureTimeFigures(Syllable).Figure(i));
           set(gca,'FontSize',14);
           hold on;
           plot([MotifDetails(Syllable).syllable.ffreq{DirectedIndices,10}],[MotifDetails(Syllable).syllable.ffreq{DirectedIndices,(i+1)}],'r.');
           plot([MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,10}],[MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,(i+1)}],'b.');
           YMax = max(max([MotifDetails(Syllable).syllable.ffreq{DirectedIndices,(i+1)}]),max([MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,(i+1)}]));
           YMin = min(min([MotifDetails(Syllable).syllable.ffreq{DirectedIndices,(i+1)}]),min([MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,(i+1)}]));

           YRange = YMax - YMin;
           YMax = YMax + 0.1*YRange;
           YMin = YMin - 0.1*YRange;
           YRange = YMax - YMin;

           XMax = max(max([MotifDetails(Syllable).syllable.ffreq{DirectedIndices,10}]),max([MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,10}]));
           XMin = min(min([MotifDetails(Syllable).syllable.ffreq{DirectedIndices,10}]),min([MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,10}]));       

           XRange = XMax - XMin;
           XMax = XMax + 0.1*XRange;
           XMin = XMin - 0.1*XRange;
           XRange = XMax - XMin;       

           axis([XMin XMax YMin YMax]);

           if (i > 2)
               set(gca,'xtick',[]);
           else
               xlabel('Time of Day','FontSize',14);
           end
           
           axis auto;
           ReScaledAxis = ReScaleAxis(axis);
           axis(ReScaledAxis);
           TitleXCoord = ReScaledAxis(1) + (ReScaledAxis(2) - ReScaledAxis(1))/2;
           TitleYCoord = ReScaledAxis(4);
           title(PlotTitles{i},'FontSize',14,'Position',[TitleXCoord TitleYCoord]);
           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot features as a function of the bout no that the motif is part of

           axes(BoutNoFigures(Syllable).Figure(i));
           set(gca,'FontSize',14);
           hold on;
           plot(([MotifDetails(Syllable).syllable.ffreq{DirectedIndices,13}] - 0.1),[MotifDetails(Syllable).syllable.ffreq{DirectedIndices,(i+1)}],'r.');
           plot(([MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,13}] + 0.1),[MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,(i+1)}],'b.');

           YMax = max(max([MotifDetails(Syllable).syllable.ffreq{DirectedIndices,(i+1)}]),max([MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,(i+1)}]));
           YMin = min(min([MotifDetails(Syllable).syllable.ffreq{DirectedIndices,(i+1)}]),min([MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,(i+1)}]));

           YRange = YMax - YMin;
           YMax = YMax + 0.1*YRange;
           YMin = YMin - 0.1*YRange;
           YRange = YMax - YMin;

           XMax = max(max([MotifDetails(Syllable).syllable.ffreq{DirectedIndices,13}]),max([MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,13}])) + 0.1;
           XMin = min(min([MotifDetails(Syllable).syllable.ffreq{DirectedIndices,13}]),min([MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,13}])) - 0.1;       

           XRange = XMax - XMin;
           XMax = XMax + 0.1*XRange;
           XMin = XMin - 0.1*XRange;
           XRange = XMax - XMin;       

           axis([XMin XMax YMin YMax]);
           
           if (i > 2)
               set(gca,'xtick',[]);
           else
               xlabel('Bout No','FontSize',14);
               for BoutNos = 1:max([MotifDetails(Syllable).syllable.ffreq{:,13}]),
                   DirectedBouts = find([MotifDetails(Syllable).syllable.ffreq{DirectedIndices,13}] == BoutNos);
                   UnDirectedBouts = find([MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,13}] == BoutNos);
                   text((BoutNos - 0.25),(YMax),['(',num2str(length(DirectedBouts)),')']);
                   text((BoutNos + 0.1),(YMax),['(',num2str(length(UnDirectedBouts)),')']);
               end
           end

           axis auto;
           ReScaledAxis = ReScaleAxis(axis);
           axis(ReScaledAxis);
           TitleXCoord = ReScaledAxis(1) + (ReScaledAxis(2) - ReScaledAxis(1))/2;
           TitleYCoord = ReScaledAxis(4);
           title(PlotTitles{i},'FontSize',14,'Position',[TitleXCoord TitleYCoord]);
           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the features as a function of the motif no. that the motif is part
    % of within the bout

           axes(MotifNoFigures(Syllable).Figure(i));
           set(gca,'FontSize',14);
           hold on;
           plot(([MotifDetails(Syllable).syllable.ffreq{DirectedIndices,14}] - 0.1),[MotifDetails(Syllable).syllable.ffreq{DirectedIndices,(i+1)}],'r.');
           plot(([MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,14}] + 0.1),[MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,(i+1)}],'b.');

           YMax = max(max([MotifDetails(Syllable).syllable.ffreq{DirectedIndices,(i+1)}]),max([MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,(i+1)}]));
           YMin = min(min([MotifDetails(Syllable).syllable.ffreq{DirectedIndices,(i+1)}]),min([MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,(i+1)}]));

           YRange = YMax - YMin;
           YMax = YMax + 0.1*YRange;
           YMin = YMin - 0.1*YRange;
           YRange = YMax - YMin;

           XMax = max(max([MotifDetails(Syllable).syllable.ffreq{DirectedIndices,14}]),max([MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,14}])) + 0.1;
           XMin = min(min([MotifDetails(Syllable).syllable.ffreq{DirectedIndices,14}]),min([MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,14}])) - 0.1;       

           XRange = XMax - XMin;
           XMax = XMax + 0.1*XRange;
           XMin = XMin - 0.1*XRange;
           XRange = XMax - XMin;       

           axis([XMin XMax YMin YMax]);
           
           if (i > 2)
               set(gca,'xtick',[]);
           else
               xlabel('Motif No','FontSize',14);
               for MotifNos = 1:max([MotifDetails(Syllable).syllable.ffreq{:,14}]),
                   DirectedBouts = find([MotifDetails(Syllable).syllable.ffreq{DirectedIndices,14}] == MotifNos);
                   UnDirectedBouts = find([MotifDetails(Syllable).syllable.ffreq{UnDirectedIndices,14}] == MotifNos);
                   text((MotifNos - 0.25),(YMax),['(',num2str(length(DirectedBouts)),')']);
                   text((MotifNos + 0.1),(YMax),['(',num2str(length(UnDirectedBouts)),')']);
               end
           end
           
           axis auto;
           ReScaledAxis = ReScaleAxis(axis);
           axis(ReScaledAxis);
           TitleXCoord = ReScaledAxis(1) + (ReScaledAxis(2) - ReScaledAxis(1))/2;
           TitleYCoord = ReScaledAxis(4);
           title(PlotTitles{i},'FontSize',14,'Position',[TitleXCoord TitleYCoord]);           
        end
    end
end
disp('Finished');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [RescaledAxis] = ReScaleAxis(OrigAxis);

XMin = OrigAxis(1);
XMax = OrigAxis(2);
YMin = OrigAxis(3);
YMax = OrigAxis(4);

XRange = XMax - XMin;
YRange = YMax - YMin;

XMin = XMin - 0.1*XRange;
XMax = XMax + 0.1*XRange;
YMin = YMin - 0.1*YRange;
YMax = YMax + 0.1*YRange;

RescaledAxis = [XMin XMax YMin YMax];