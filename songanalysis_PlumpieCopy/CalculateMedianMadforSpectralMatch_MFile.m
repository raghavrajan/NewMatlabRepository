function [WinMedian, WinMAD] = CalculateMedianMadforSpectralMatch_MFile(TempS, WMotif)

WinMedian = zeros((size(TempS,2) - size(WMotif,2) + 1), 1);
WinMAD = zeros((size(TempS,2) - size(WMotif,2) + 1), 1);
                
for ColNo = 1:(size(TempS,2) - size(WMotif, 2) + 1),
    StartIndex = ((ColNo - 1)*size(WMotif,1)) + 1;
    WinIndices = StartIndex:1:(StartIndex + size(WMotif,1)*size(WMotif,2) - 1);
    WinMedian(ColNo) = median(TempS(WinIndices));
    WinMAD(ColNo) = mad(TempS(WinIndices));
end
