function [Dir_p, UnDir_p] = MA_PlotDirVsUnDir(DirMedians, UnDirMedians, YLabelString, OutputDir, FileNameString, TitleString, varargin) 

if (nargin > 6)
    PlotYesNo = varargin{1};
else
    PlotYesNo = 1;
end

if (PlotYesNo == 1)
    BarPlotFigure = figure;
    set(gcf, 'Color', 'w');
    hold on;
    Colors = ['rb'];
    for i = 1:size(DirMedians, 2),
        DirBar(i) = bar(i, mean(DirMedians(:,i)));
        set(DirBar(i), 'FaceColor', 'none', 'EdgeColor', Colors(i));
        errorbar(i, mean(DirMedians(:,i)), std(DirMedians(:,i))/sqrt(size(DirMedians, 1)), [Colors(i), 'o']);
    end
    plot(repmat([1.2 1.8], size(DirMedians,1), 1)', (DirMedians)', 'ko-')

    for i = 1:size(UnDirMedians, 2),
        UnDirBar(i) = bar(i+3, mean(UnDirMedians(:,i)));
        set(UnDirBar(i), 'FaceColor', 'none', 'EdgeColor', Colors(i));
        errorbar(i+3, mean(UnDirMedians(:,i)), std(UnDirMedians(:,i))/sqrt(size(UnDirMedians, 1)), [Colors(i), 'o']);
    end
    plot(repmat([4.2 4.8], size(UnDirMedians,1), 1)', (UnDirMedians)', 'ko-');
end

Dir_p = signrank(DirMedians(:,1), DirMedians(:,2));
UnDir_p = signrank(UnDirMedians(:,1), UnDirMedians(:,2));

if (PlotYesNo == 1)
    PlotAxis = [0.5 5.5 min([0 (0.95*min([UnDirMedians(:); DirMedians(:)]))]) max([0 (1.05*max([UnDirMedians(:); DirMedians(:)]))])];
    SetLabelsForFigure(BarPlotFigure, 'YLabel', YLabelString, 'XTickNums', '1 2 4 5', 'XTickLabel', 'Dir Undir Dir Undir', 'Axis', num2str(PlotAxis));
    SetLabelsForFigure(BarPlotFigure, 'Text', ['Pre (p=', num2str(Dir_p), ')'], ['0.9 ', num2str(PlotAxis(4)/1.01)], 'Text', ['Post (p=', num2str(UnDir_p), ')'], ['3.9 ', num2str(PlotAxis(4)/1.01)]);
    saveas(BarPlotFigure, [OutputDir, TitleString, '.', FileNameString, '.DirVsUnDir.BarPlot.png'], 'png');
end