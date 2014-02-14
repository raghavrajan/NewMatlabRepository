function [] = PlotRaster(Raster, varargin)

% =========================================================================
% PlotRaster function to plot raster - RR
% Usage: 
%      1) PlotRaster(Raster)
%          Raster is a two column array where the first column corresponds
%          to spike times and the second column corresponds to the trial
%          index
%
%      2) PlotRaster(Raster, Colour)
%          Colour refers to the colour for the ticks
%
%      3) PlotRaster(Raster, Colour, LineWidth)
%          LineWidth is a scalar that is used for the width of each tick
%
%      4) PlotRaster(Raster, Colour, LineWidth, RasterIncrement)
%          RasterIncrement is a scalar that is added to the trial indices
%          to offset the entire raster. (Useful for plotting multiple
%          rasters from different Raster arrays on the same plot)
%
%      5) PlotRaster(Raster, Colour, LineWidth, RasterIncrement, MaxTrials)
%          MaxTrials is a scalar that can be used to plot only a subset of
%          the trials for a given raster
%
%==========================================================================

if (nargin > 1)
    Colour = varargin{1};
end
if (nargin > 2)
    LW = varargin{2};
end
if (nargin > 3)
    RasterIncrement = varargin{3};
end
if (nargin > 4)
    MaxRaster = varargin{4};
end

if (exist('MaxRaster', 'var'))
    Raster = Raster(find(Raster(:,2) <= MaxRaster),:);
end

if (exist('RasterIncrement', 'var'))
    Raster(:,2) = Raster(:,2) + RasterIncrement;
end

hold on;
Raster = line([Raster(:,1)'; Raster(:,1)'], [(Raster(:,2)' - 0.25); (Raster(:,2)' + 0.25)], 'Color', 'k', 'LineWidth', 0.25);

if (exist('Colour', 'var'))
    set(Raster, 'Color', Colour);
end

if (exist('LW', 'var'))
    set(Raster, 'LineWidth', LW);
end
    
%MarkerString = repmat('|',size(Raster,1),1);
%text(Raster(:,1),Raster(:,2),MarkerString,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',2,'FontName','FixedWidth','FontUnits','pixels','Color', 'k','FontWeight','bold');