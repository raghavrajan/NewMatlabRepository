function [BeeSwarmHandle] = BeeSwarmPlot(XData, YData, AxesHandle, varargin)

% For marker size of 2, I can put points them at about 0.07 apart. I can
% also spread them out by 0.28 on either side of the category point. So I
% should decide how many points there are and then spread them out. Given
% the values I have put above, I can manage about 4 points on either side
% of the center point. So, I should decide on whether there are more than 9
% points at a given location. If there are, then I make the marker as .,
% otherwise I make the marker as o.

if (nargin > 3)
    MarkerColourValue = varargin{1};
else
    MarkerColourValue = 'k';
end

if (nargin > 4)
    FilledOrNot = varargin{2};
else
    FilledOrNot = 0;
end

if (FilledOrNot == 1)
    MarkerFaceColourValue = MarkerColourValue;
else
    MarkerFaceColourValue = 'none';
end
axes(AxesHandle);
hold on;

% So, first check whether the mode of the data is more than 9 or not.
DataMode = mode(YData);

for i = 1:length(DataMode),
    NumModalValues(i) = length(find(YData == DataMode(i)));
end

if (max(NumModalValues) > 9)
    MarkerSymbol = '.';
    switch (max(NumModalValues))
        case {10 11 12 13 14 15}
            MarkerSizeValue = 6;
            InterMarkerDistance = 0.04;
            
        otherwise
            MarkerSizeValue = 4;
            InterMarkerDistance = 0.04;
    end
else
    MarkerSymbol = 'o';
    switch (max(NumModalValues))
        case {3 2 1}
            MarkerSizeValue = 6;
            InterMarkerDistance = 0.15;
            
        case {7 6 5 4}
            MarkerSizeValue = 4;
            InterMarkerDistance = 0.1;
            
        otherwise
            MarkerSizeValue = 2;
            InterMarkerDistance = 0.07;
    end
end

% Now spread out the values accordingly
UniqueYValues = unique(YData);
for i = UniqueYValues(:)',
    Matches = find(YData == i);
    if (mod(length(Matches),2) == 1) % odd - then one point in the center and the rest spread out plus and minus of the center
        XSpreadData = XData + ((-(((length(Matches)-1)/2) * InterMarkerDistance)):InterMarkerDistance:(((length(Matches) - 1)/2) * InterMarkerDistance));
    else
        % Now it is an even number of points. So then I have to spread out
        % half the points on one side of the center and half on the other
        % side to make it even
        XSpreadData = XData + ((-((length(Matches)/2) * InterMarkerDistance) + InterMarkerDistance/2):InterMarkerDistance:(((length(Matches)/2) * InterMarkerDistance) - InterMarkerDistance/2));
    end
    plot(XSpreadData, YData(Matches), ['k', MarkerSymbol], 'MarkerSize', MarkerSizeValue, 'MarkerFaceColor', MarkerColourValue, 'Color', MarkerColourValue, 'MarkerFaceColor', MarkerFaceColourValue);
end





