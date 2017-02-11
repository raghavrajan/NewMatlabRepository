function [] = CompareFiringRateEstimates(SpikeTrain, FileLength)

% =========================================================================
% CompareFiringRateEstimates is a script to compare different ways of 
% calculating firing rate as a continuous function of time
%
% Usage: CompareFiringRateEstimates(SpikeTrain, FileLength)
%
% Inputs - 1) SpikeTrain is a cell array with each cell having one trial of
%             spike times
%          2) FileLength has two values indicating the start and end of the 
%             period for which spike times are stored in SpikeTrain. It is
%             assumed that all the entries in SpikeTrain are for periods of
%             equal length.
%
% 
% Raghav - 20th September 2011
% =========================================================================

% First, using the instantaneous firing rate which is defined as
%
% F(t) = 1/(t(i) - t(i-1)) over the interval t(i-1)<t<=t(i)
%
% where F(t) is the firing rate, t(i) and t(i-1) are the times of the ith
% and (i-1)th spike respectively

IFRdt = 0.001; % in seconds

for j = 1:length(SpikeTrain),
    Time = FileLength(1):IFRdt:FileLength(2);
    IFR(j,:) = CalculateIFR(SpikeTrain{j} + abs(FileLength(1)), Time);
end

% Now use a bin size of 10ms and calculate PSTH
Edges_10ms = FileLength(1):0.01:FileLength(2);
for j = 1:length(SpikeTrain),
    PST_10ms(j,:) = histc(SpikeTrain{j}, Edges_10ms);
end
PST_10ms = PST_10ms/0.01;

% Now use a bin size of 5ms and calculate PSTH
Edges_5ms = FileLength(1):0.005:FileLength(2);
for j = 1:length(SpikeTrain),
    PST_5ms(j,:) = histc(SpikeTrain{j}, Edges_5ms);
end
PST_5ms = PST_5ms/0.005;

% Now use Bars estimate
Fit1 = barsP(mean(PST_5ms), [Edges_5ms(1) Edges_5ms(end)], 372);

% Now use Gaussian smoothing with 5ms window
GaussianLen = 4;
Width = 0.005;
Fs = 1/(Edges_5ms(2) - Edges_5ms(1));

XGauss = 1:1:(1 + round(2 * GaussianLen * Width * (Fs)));
XGauss = XGauss - (length(XGauss) + 1)/2;
GaussWin = (1/((Width * Fs) * sqrt(2 * pi))) * exp(-(XGauss.*XGauss)/(2 * (Width * Fs) * (Width * Fs)));

dt = 1/Fs;

figure;
RasterAxis = subplot(2,1,1);
for i = 1:length(SpikeTrain),
    PlotRaster([SpikeTrain{i} ones(length(SpikeTrain{i}), 1)*i], 'k', 1);
end

PSTAxis = subplot(2,1,2);
hold on;
plot(Time, mean(IFR), 'r');
plot(Edges_10ms, mean(PST_10ms), 'b');
plot(Edges_5ms, mean(PST_5ms), 'c');
plot(Edges_5ms, Fit1.mean, 'm');
legend('IFR', '1ms', '10ms', '5ms', '5ms BARS');
%PST_1msBar = bar(Edges_1ms, mean(PST_1ms));
%set(PST_1msBar, 'FaceColor', 'none', 'EdgeColor', 'k', 'BarWidth', 0.2);
disp('Finished all plots');