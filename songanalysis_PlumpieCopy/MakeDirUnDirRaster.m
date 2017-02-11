function [] = MakeDirUnDirRaster()

MainFigure = figure;
set(MainFigure, 'Color', 'w');

PSTPlot = axes('Position',[0.15 0.05 0.55 0.1]);
SpontActivityPlot = axes('Position',[0.15 0.18 0.55 0.1]);
UnDirRasterPlot = axes('Position',[0.15 0.31 0.55 0.2]);
DirRasterPlot = axes('Position',[0.15 0.54 0.55 0.2]);
UnDirMotifPlot = axes('Position',[0.15 0.77 0.55 0.08]);
DirMotifPlot = axes('Position',[0.15 0.88 0.55 0.08]);

SpontISIPlot = axes('Position', [0.8 0.2 0.15 0.12]);
UnDirISIPlot = axes('Position', [0.8 0.41 0.15 0.12]);
DirISIPlot = axes('Position', [0.8 0.64 0.15 0.12]);
WaveformPlot = axes('Position', [0.8 0.8 0.15 0.1]);

