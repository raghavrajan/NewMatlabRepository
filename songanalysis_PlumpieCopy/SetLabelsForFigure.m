function [] = SetLabelsForFigure(Fig, varargin)

% Script to scale the axis and put xlabel, ylabels, set the position, etc.


%=========================================================================
% Defaults
figure(Fig);
set(Fig, 'Color', 'w');
set(Fig, 'PaperPositionMode', 'auto');
SubPlotYesNo = 0;

if (nargin > 1)
    
    % Find out whether there are subplots or not
    Index = find(cellfun(@length, strfind(varargin, 'Subplot')));
    if (~isempty(Index))
        SubplotSize = cell2mat(textscan(varargin{Index + 1}, '%n%n'));
        SubPlotYesNo = 1;
        for i = 1:SubplotSize(1)*SubplotSize(2),
            subplot(SubplotSize(1), SubplotSize(2), i);
            set(gca, 'FontSize', 12, 'FontWeight', 'bold');
            set(gca, 'Box', 'off');
        end
    else
        set(gca, 'FontSize', 12, 'FontWeight', 'bold');
        set(gca, 'Box', 'off');
    end
    
    % Set axis
    Index = find(cellfun(@length, strfind(varargin, 'Axis')));
    if (~isempty(Index))
        PlotAxis = cell2mat(textscan(varargin{Index + 1}, '%n%n%n%n'));
        if (SubPlotYesNo == 0)
            axis(PlotAxis);
        else
            for i = 1:SubplotSize(1)*SubplotSize(2),
                subplot(SubplotSize(1), SubplotSize(2), i);
                axis(PlotAxis)
            end
        end
    end    
    
    % Set X-axis ticks and their labels
    Index = find(cellfun(@length, strfind(varargin, 'XTickNums')));
    if (~isempty(Index))
        XTickNums = cell2mat(textscan(varargin{Index + 1}, '%n', 'DeLimiter', ' '));
        Index = find(cellfun(@length, strfind(varargin, 'XTickLabel')));
        if (~isempty(Index))
            XTickLabels = textscan(varargin{Index + 1}, '%s', 'DeLimiter', ' ');
            XTickLabels = XTickLabels{1};
            if (SubPlotYesNo == 0)
                set(gca, 'XTick', XTickNums, 'XTickLabel', XTickLabels);
            else
                for i = 1:SubplotSize(1)*SubplotSize(2),
                    subplot(SubplotSize(1), SubplotSize(2), i);
                    set(gca, 'XTick', XTickNums, 'XTickLabel', XTickLabels);
                end
            end
        end
    end
    
    % Set xlabel
    Index = find(cellfun(@length, strfind(varargin, 'XLabel')));
    if (~isempty(Index))
        if (SubPlotYesNo == 0)
            xlabel(varargin{Index + 1}, 'FontSize', 14, 'FontWeight', 'bold');
        else
            if (SubplotSize(2) > 1)
                LastRowPlots = (SubplotSize(1)*(SubplotSize(2) - 1) + 1):1:(SubplotSize(1)*SubplotSize(2)); 
            else
                LastRowPlots = SubplotSize(1)*SubplotSize(2);
            end
            subplot(SubplotSize(1), SubplotSize(2), median(floor(LastRowPlots)));
            xlabel(varargin{Index + 1}, 'FontSize', 14, 'FontWeight', 'bold');
        end
    end
    
    % Set ylabel
    Index = find(cellfun(@length, strfind(varargin, 'YLabel')));
    if (~isempty(Index))
        if (SubPlotYesNo == 0)
            ylabel(varargin{Index + 1}, 'FontSize', 14, 'FontWeight', 'bold');
        else
            FirstColPlots = 1:SubplotSize(2):(SubplotSize(1)*SubplotSize(2));
            subplot(SubplotSize(1), SubplotSize(2), median(floor(FirstColPlots)));
            ylabel(varargin{Index + 1}, 'FontSize', 14, 'FontWeight', 'bold');
        end
    end
    
    % set titles
    Index = find(cellfun(@length, strfind(varargin, 'Title')));
    if (~isempty(Index))
        Titles = textscan(varargin{Index + 1}, '%s', 'DeLimiter', ' ');
        Titles = Titles{1};
        if (SubPlotYesNo == 0)
            title(Titles, 'FontSize', 14, 'FontWeight', 'bold');
        else
            for j = 1:SubplotSize(1)*SubplotSize(2),
                subplot(SubplotSize(1), SubplotSize(2), j);
                title(Titles{j}, 'FontSize', 14, 'FontWeight', 'bold');
            end
        end
    end
    
    % set Legend
    Index = find(cellfun(@length, strfind(varargin, 'Legend')));
    if (~isempty(Index))
        Legends = textscan(varargin{Index + 1}, '%s', 'DeLimiter', ' ');
        Legends = Legends{1};
        if (SubPlotYesNo == 0)
            legend(Legends);
            legend('boxoff');
            legend('Location', varargin{Index + 2});
        else
            for j = 1:SubplotSize(1)*SubplotSize(2),
                subplot(SubplotSize(1), SubplotSize(2), j);
                legend(Legends);
                legend('boxoff');
                legend('Location', varargin{Index + 2});
            end
        end
    end
    
    % Write text
    Index = find(cellfun(@length, strfind(varargin, 'Text')));
    if (~isempty(Index))
        if (SubPlotYesNo == 0)
            for i = 1:length(Index),
                String = varargin{Index(i) + 1};
                TextCoords = cell2mat(textscan(varargin{Index(i) + 2}, '%n%n'));
                text(TextCoords(1), TextCoords(2), String, 'FontSize', 12, 'FontWeight', 'bold');
            end
        else
            for j = 1:SubplotSize(1)*SubplotSize(2),
                subplot(SubplotSize(1), SubplotSize(2), j);
                for i = 1:length(Index),
                    String = varargin{Index(i) + 1};
                    TextCoords = cell2mat(textscan(varargin{Index(i) + 2}, '%n%n'));
                    text(TextCoords(1), TextCoords(2), String, 'FontSize', 12, 'FontWeigth', 'bold');
                end
            end
        end
    end
    
    % Draw a line parallel to x-axis
    Index = find(cellfun(@length, strfind(varargin, 'LineX')));
    if (~isempty(Index))
        if (SubPlotYesNo == 0)
            for i = 1:length(Index),
                Coords = cell2mat(textscan(varargin{Index(i) + 1}, '%n'));
                Temp = axis;
                plot(Temp(1:2), ones(2,1)*Coords, 'k--', 'LineWidth', 2);
            end
        else
            for j = 1:SubplotSize(1)*SubplotSize(2),
                subplot(SubplotSize(1), SubplotSize(2), j);
                for i = 1:length(Index),
                    Coords = cell2mat(textscan(varargin{Index(i) + 1}, '%n'));
                    Temp = axis;
                    plot(Temp(1:2), ones(2,1)*Coords, 'k--', 'LineWidth', 2);
                end
            end
        end
    end
    
    % Draw a line parallel to y-axis
    Index = find(cellfun(@length, strfind(varargin, 'LineY')));
    if (~isempty(Index))
        if (SubPlotYesNo == 0)
            for i = 1:length(Index),
                Coords = cell2mat(textscan(varargin{Index(i) + 1}, '%n'));
                Temp = axis;
                plot(ones(2,1)*Coords, Temp(3:4), 'k--', 'LineWidth', 2);
            end
        else
            for j = 1:SubplotSize(1)*SubplotSize(2),
                subplot(SubplotSize(1), SubplotSize(2), j);
                for i = 1:length(Index),
                    Coords = cell2mat(textscan(varargin{Index(i) + 1}, '%n'));
                    Temp = axis;
                    plot(ones(2,1)*Coords, Temp(3:4), 'k--', 'LineWidth', 2);
                end
            end
        end
    end

end