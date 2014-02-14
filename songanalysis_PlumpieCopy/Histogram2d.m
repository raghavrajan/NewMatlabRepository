function [Values] = Histogram2d(X, XAxis, Y, YAxis, XBins, YBins)

XBinEdges = linspace(XAxis(1), XAxis(2), (XBins + 1));
YBinEdges = linspace(YAxis(1), YAxis(2), (YBins + 1));

Values = zeros(XBins, YBins);

for i = 2:length(XBinEdges),
    Indices = find((X >= XBinEdges(i-1)) & (X < XBinEdges(i)));
    for j = 2:length(YBinEdges),
        Values((i-1),(j-1)) = length(find((Y(Indices) >= YBinEdges(j-1)) & (Y(Indices) < YBinEdges(j))));
    end
end
    
disp('Done');