function [] = PlotISI(SpikeRaster, MainFigure, ISIPlot, TitleString)

figure(MainFigure);
axes(ISIPlot);

ISI = diff(SpikeRaster(:,1));

ISI(find(ISI < 0)) = 2;

ISI = ISI * 1000;

Edges = 0:1:1000;

ISIs = histc(ISI, Edges);
bar(Edges, ISIs, 'histc');
axis tight;
temp = axis;
set(gca, 'XScale', 'log');
axis([1 1000 0 (temp(4) * 1.2)]);
set(gca, 'Box', 'off');
set(gca, 'XTick', [1 10 100 1000], 'XTickLabel', [1 10 100 1000]);
title(TitleString);
xlabel('ISI ms');
ylabel('Count');
