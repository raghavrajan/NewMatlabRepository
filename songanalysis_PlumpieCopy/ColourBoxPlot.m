function [ColourBoxPlotHandle] = ColourBoxPlot(XData, YData, LineColour, Width, FilledOrNot, varargin)

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

if (isempty(YData))
    disp('No Y data');
    ColourBoxPlotHandle = [];
    return;
end


MinNumberPoints = 3;
YData = sort(YData);

hold on;
WidthIncreaseFactor = 1;
WhiskerIncreaseFactor = 0.75;
ActualLineWidth = 1.5;

%% First plot median
ColourBoxPlotHandle(1) = plot(XData + [-Width/2 Width/2], ones(1,2) * median(YData), LineColour, 'LineWidth', ActualLineWidth);

%% Next plot the 25th and 75th percentile - make the lines a little longer to make a trapezium
if (length(YData) >= MinNumberPoints)
    ColourPlotHandle(2) = patch([(ones(1,2) * (XData + WidthIncreaseFactor*(-Width/2))) (ones(1,2) * (XData + WidthIncreaseFactor*Width/2))], [prctile(YData, [25 75]) prctile(YData, [75 25])], LineColour, 'LineWidth', ActualLineWidth);
    switch FilledOrNot
        case 'filled'
            set(ColourPlotHandle(2), 'EdgeColor', LineColour, 'FaceAlpha', 0.1);
        otherwise
            set(ColourPlotHandle(2), 'EdgeColor', LineColour, 'FaceColor', 'none');
    end
else
    plot(ones(size(YData))*XData, YData, [LineColour, 'o'], 'MarkerFaceColor', LineColour, 'LineWidth', ActualLineWidth);
end

%% Now plot whiskers - largest/smallests points within +/- 1.5 * IQR
if (length(YData) >= MinNumberPoints)
    Thresholds = [(prctile(YData, 25) - 1.5*iqr(YData)) (prctile(YData, 75) + 1.5*iqr(YData))];
    
    UpperWhiskerExtreme = YData(find((YData > prctile(YData, 75)) & (YData <= Thresholds(2))));
    if (~isempty(UpperWhiskerExtreme))
        UpperWhiskerExtreme = UpperWhiskerExtreme(end);
        ColourBoxPlotHandle(end+1) = plot(XData*ones(1,2), [prctile(YData, 75) UpperWhiskerExtreme], LineColour, 'LineWidth', ActualLineWidth);
        ColourBoxPlotHandle(end+1) = plot(XData + (WhiskerIncreaseFactor*[-Width/2 Width/2]), ones(1,2) * UpperWhiskerExtreme, LineColour, 'LineWidth', ActualLineWidth);
    end

    LowerWhiskerExtreme = YData(find((YData >= Thresholds(1)) & (YData < prctile(YData, 25))));
    if (~isempty(LowerWhiskerExtreme))
        LowerWhiskerExtreme = LowerWhiskerExtreme(1);
        ColourBoxPlotHandle(end+1) = plot(XData*ones(1,2), [prctile(YData, 25) LowerWhiskerExtreme], LineColour, 'LineWidth', ActualLineWidth);
        ColourBoxPlotHandle(end+1) = plot(XData + (WhiskerIncreaseFactor*[-Width/2 Width/2]), ones(1,2) * LowerWhiskerExtreme, LineColour, 'LineWidth', ActualLineWidth);
    end
end

if (PlotOutliersOrNot == 1)
    OutlierLineColour = 'r';
    %% Now plot the outliers as crosses (filled) - outside of +/- 1*IQR from 25% and 75%
    if (length(YData) >= MinNumberPoints)
        Outliers = find((YData < Thresholds(1)) | (YData > Thresholds(2)));
        if (~isempty(Outliers))
            ColourBoxPlotHandle(end+1) = plot(XData * ones(size(Outliers)), YData(Outliers), [OutlierLineColour, '+'], 'MarkerSize', 6, 'MarkerFaceColor', LineColour, 'LineWidth', ActualLineWidth);
        end
    end
end

%% Plot raw data points too if specified in the input
if (PlotRawDataOrNot == 1)
    plot(ones(size(YData))*XData + (Width*rand(size(YData)) - Width/2), YData, 'ko', 'MarkerSize', 1, 'MarkerFaceColor', [0.5 0.5 0.5], 'Color', [0.5 0.5 0.5]);
end

%% Now also write the # of data points above the upper whisker
if (length(YData) >= MinNumberPoints)
    if (isempty(UpperWhiskerExtreme))
        UpperWhiskerExtreme = prctile(YData, 75);
    end
else
    UpperWhiskerExtreme = max(YData);
end
% text(XData, 1.1 * UpperWhiskerExtreme, ['(', num2str(length(YData)), ')'], 'FontSize', 16, 'HorizontalAlignment', 'center');
