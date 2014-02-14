function [OptimalBinSize] = EstimateOptimalBinSize(SpikeTrain, FileLength)

% =========================================================================
% EstimateOptimalBinSize is a script to estimate the optimal bin size for a
% set of data using the method suggested by Shimazaki and Shinomoto 2007
%
% Usage: OptimalBinSize = EstimateOptimalBinSize(SpikeTrain, FileLength)
%
% Inputs - 1) SpikeTrain is a cell array with each cell having one trial of
%             spike times
%          2) FileLength has two values indicating the start and end of the 
%             period for which spike times are stored in SpikeTrain. It is
%             assumed that all the entries in SpikeTrain are for periods of
%             equal length.
% Outputs - OptimalBinSize in sec
% 
% Raghav - 20th September 2011
% =========================================================================

SpikeTimes = [];
for j = 1:length(SpikeTrain),
    SpikeTimes = [SpikeTimes; SpikeTrain{j}];
end

Index = 0;
for NoofBins = 2:2000,
    Index = Index + 1;
    Edges = linspace(FileLength(1), FileLength(2), NoofBins);
    PST = histc(SpikeTimes, Edges);
    PST = PST(1:end-1);
    k = sum(PST)/NoofBins;
    v = sum((PST - k).*(PST - k))/NoofBins;
    BinWidth = Edges(2) - Edges(1);
    C(Index,:) = [BinWidth (2*k - v)/(BinWidth * BinWidth * NoofBins)];
end
disp('Calculated optimal bin size');