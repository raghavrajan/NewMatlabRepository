function [ColourBoxPlotHandle] = ColourBoxPlot_XAxis(XData, YData, LineColour, Width, FilledOrNot, varargin)

% =========================================================================
% Script for creating a box plot for specific x and y data - the boxplot
% that is in-built in MATLAB only plots for a certain set of data and the
% associated group. I want to be able to specify one data set and do a box
% plot of it and then do another set if necessary etc. separately.
% Boxplot will have median, 25%, 75% and whiskers at 1*IQR. Points outside
% of that should be plotted individually.
% This script will also decide on whether to plot as a box-plot (if there
% are 3 or more points), otherwise just plot the individual points and the
% median
% =========================================================================

if (nargin > 5)
    PlotOutliersOrNot = varargin{1};
else
    PlotOutliersOrNot = 0;
end

if (nargin > 6)
    PlotRawDataOrNot = varargin{2};
else
    PlotRawDataOrNot = 0;
end

if (isempty(XData))
    disp('No X data');
    ColourBoxPlotHandle = [];
    return;
end


MinNumberPoints = 3;
XData = sort(XData);

hold on;
WidthIncreaseFactor = 1;
WhiskerIncreaseFactor = 0.75;
ActualLineWidth = 1.5;

%% First plot median
ColourBoxPlotHandle(1) = plot(ones(1,2) * median(XData), YData + [-Width/2 Width/2], LineColour, 'LineWidth', ActualLineWidth);

%% Next plot the 25th and 75th percentile - make the lines a little longer to make a trapezium
if (length(XData) >= MinNumberPoints)
    ColourPlotHandle(2) = patch([prctile(XData, [25 75]) prctile(XData, [75 25])], [(ones(1,2) * (YData + WidthIncreaseFactor*(-Width/2))) (ones(1,2) * (YData + WidthIncreaseFactor*Width/2))], LineColour, 'LineWidth', ActualLineWidth);
    switch FilledOrNot
        case 'filled'
            set(ColourPlotHandle(2), 'EdgeColor', LineColour, 'FaceAlpha', 0.1);
        otherwise
            set(ColourPlotHandle(2), 'EdgeColor', LineColour, 'FaceColor', 'none');
    end
else
    plot(XData, ones(size(XData))*YData, [LineColour, 'o'], 'MarkerFaceColor', LineColour, 'LineWidth', ActualLineWidth);
end

% %% Now plot whiskers - largest/smallests points within +/- 1.5 * IQR
if (length(XData) >= MinNumberPoints)
    Thresholds = [(prctile(XData, 25) - 1.5*iqr(XData)) (prctile(XData, 75) + 1.5*iqr(XData))];
    
    UpperWhiskerExtreme = XData(find((XData > prctile(XData, 75)) & (XData <= Thresholds(2))));
    if (~isempty(UpperWhiskerExtreme))
        UpperWhiskerExtreme = UpperWhiskerExtreme(end);
        ColourBoxPlotHandle(end+1) = plot([prctile(XData, 75) UpperWhiskerExtreme], YData*ones(1,2), LineColour, 'LineWidth', ActualLineWidth);
        ColourBoxPlotHandle(end+1) = plot(ones(1,2) * UpperWhiskerExtreme, YData + (WhiskerIncreaseFactor*[-Width/2 Width/2]), LineColour, 'LineWidth', ActualLineWidth);
    end

    LowerWhiskerExtreme = XData(find((XData >= Thresholds(1)) & (XData < prctile(XData, 25))));
    if (~isempty(LowerWhiskerExtreme))
        LowerWhiskerExtreme = LowerWhiskerExtreme(1);
        ColourBoxPlotHandle(end+1) = plot([prctile(XData, 25) LowerWhiskerExtreme], YData*ones(1,2), LineColour, 'LineWidth', ActualLineWidth);
        ColourBoxPlotHandle(end+1) = plot(ones(1,2) * LowerWhiskerExtreme, YData + (WhiskerIncreaseFactor*[-Width/2 Width/2]), LineColour, 'LineWidth', ActualLineWidth);
    end
end

if (PlotOutliersOrNot == 1)
    OutlierLineColour = 'r';
    %% Now plot the outliers as crosses (filled) - outside of +/- 1*IQR from 25% and 75%
    if (length(XData) >= MinNumberPoints)
        Outliers = find((XData < Thresholds(1)) | (XData > Thresholds(2)));
        if (~isempty(Outliers))
            ColourBoxPlotHandle(end+1) = plot(XData(Outliers), YData * ones(size(Outliers)), [OutlierLineColour, '+'], 'MarkerSize', 6, 'MarkerFaceColor', LineColour, 'LineWidth', ActualLineWidth);
        end
    end
end

%% Plot raw data points too if specified in the input
if (PlotRawDataOrNot == 1)
    plot(XData, ones(size(XData))*YData + (Width*rand(size(XData)) - Width/2), 'ko', 'MarkerSize', 2, 'MarkerFaceColor', [0.5 0.5 0.5], 'Color', [0.5 0.5 0.5]);
end

%% Now also write the # of data points above the upper whisker
if (length(XData) >= MinNumberPoints)
    if (isempty(UpperWhiskerExtreme))
        UpperWhiskerExtreme = prctile(XData, 75);
    end
else
    UpperWhiskerExtreme = max(XData);
end
% text(1.1 * UpperWhiskerExtreme, YData, ['(', num2str(length(XData)), ')'], 'FontSize', 16, 'HorizontalAlignment', 'center');
