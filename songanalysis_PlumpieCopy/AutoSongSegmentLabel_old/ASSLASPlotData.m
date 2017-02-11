function [] = ASSLASPlotData(handles, varargin)

if (nargin > 1)
    DisplayXLabel = varargin{1};
end

if (nargin > 2)
    Color = varargin{2};
end

if (nargin > 3)
    Symbol = varargin{3};
end

if (nargin > 4)
    Size = varargin{4};
end

if (nargin > 5)
    AxisName = varargin{5};
end

if (nargin > 6)
    ClearAxis = varargin{6};
else
    ClearAxis = 1;
end

if (nargin > 7)
    Indices = varargin{7};
else
    Indices = 1:1:size(handles.DataStruct.FeatValues, 1);
end

if (ClearAxis == 1)
    cla(AxisName);
end

axes(AxisName);
Temp = plot(handles.DataStruct.FeatValues(Indices,handles.ASSLAS.XVal), handles.DataStruct.FeatValues(Indices,handles.ASSLAS.YVal), 'k+', 'MarkerSize', 2);
set(Temp, 'Color', Color, 'Marker', Symbol, 'MarkerSize', Size);

hold on;
if (DisplayXLabel == 1)
    xlabel(handles.ASSLAS.FeatNames{handles.ASSLAS.XVal}, 'FontSize', 12, 'FontWeight', 'bold');
end
ylabel(handles.ASSLAS.FeatNames{handles.ASSLAS.YVal}, 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
axis tight;
